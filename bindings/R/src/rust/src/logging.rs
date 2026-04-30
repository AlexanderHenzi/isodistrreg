use crossbeam_channel::{Receiver, RecvTimeoutError, Sender, unbounded};
use extendr_api::print_r_output;
use indicatif::{ProgressBar, ProgressDrawTarget, ProgressStyle, TermLike};
use isodistrreg::ProgressTracker;
use std::time::Duration;

unsafe extern "C" {
    fn R_FlushConsole();
}

/// Hard-coded console width. R's `getOption("width")` defaults to 80 and most environments
/// (terminal R, RStudio, knitr) match or exceed that. Indicatif uses this both to size the bar
/// and to compute trailing padding; if the real terminal is narrower the bar wraps and leaves
/// stale lines, so we keep our rendered template comfortably under this limit.
const R_CONSOLE_WIDTH: u16 = 80;

/// Routes indicatif's draw calls through R's console output stream so the progress bar appears in
/// plain terminal R, RStudio, knitr, and anywhere `Rprintf` lands — not just on stdout.
///
/// This works because indicatif renders single-line bars by writing `\r` + content + padding +
/// flush; the `move_cursor_*` calls reduce to no-ops with `n = 0`. RStudio's console is fine
/// with `\r` line rewrites, so we don't need ANSI cursor support.
#[derive(Debug)]
struct RTerm;

impl TermLike for RTerm {
    fn width(&self) -> u16 {
        R_CONSOLE_WIDTH
    }
    fn move_cursor_up(&self, _n: usize) -> std::io::Result<()> {
        Ok(())
    }
    fn move_cursor_down(&self, _n: usize) -> std::io::Result<()> {
        Ok(())
    }
    fn move_cursor_right(&self, _n: usize) -> std::io::Result<()> {
        Ok(())
    }
    fn move_cursor_left(&self, _n: usize) -> std::io::Result<()> {
        Ok(())
    }
    fn write_line(&self, s: &str) -> std::io::Result<()> {
        print_r_output(s);
        print_r_output("\n");
        Ok(())
    }
    fn write_str(&self, s: &str) -> std::io::Result<()> {
        print_r_output(s);
        Ok(())
    }
    fn clear_line(&self) -> std::io::Result<()> {
        print_r_output("\r");
        print_r_output(" ".repeat(R_CONSOLE_WIDTH as usize));
        print_r_output("\r");
        Ok(())
    }
    fn flush(&self) -> std::io::Result<()> {
        unsafe { R_FlushConsole() };
        Ok(())
    }
}

/// Rate at which the main thread wakes up to drain pending progress messages. Matches indicatif's
/// 8 Hz draw rate, so we never sit on a redrawable update for longer than one frame.
const PUMP_TICK: Duration = Duration::from_millis(125);

/// Update from a worker thread to the main-thread renderer. Worker threads must never invoke
/// `print_r_output` directly — `RTerm` calls into R, and R's console FFI is not thread-safe.
pub enum ProgressMsg {
    SetTotal(u64),
    Increment(u64),
}

/// Worker-side handle. `set_total` and `increment` only push messages onto a channel; the
/// matching [`ProgressPump`] (running on the R-callable thread) drains the channel and forwards
/// updates to the indicatif `ProgressBar`. This keeps every R FFI call on the main thread.
pub struct IndicatifProgress {
    tx: Sender<ProgressMsg>,
}

impl ProgressTracker for IndicatifProgress {
    fn set_total(&self, n: usize) {
        let _ = self.tx.send(ProgressMsg::SetTotal(n as u64));
    }
    fn increment(&self) {
        let _ = self.tx.send(ProgressMsg::Increment(1));
    }
    /// No-op on the worker side — the bar is finalized by [`ProgressPump::finish`] on the main
    /// thread, after all pending messages have been drained.
    fn finish(&self) {}
}

/// Main-thread side of the progress channel. Owns the `ProgressBar` (which renders through
/// [`RTerm`] → `print_r_output`) and must only be driven from the thread R called us on.
pub struct ProgressPump {
    bar: ProgressBar,
    rx: Receiver<ProgressMsg>,
}

impl ProgressPump {
    pub fn new() -> (IndicatifProgress, Self) {
        let bar = ProgressBar::with_draw_target(
            Some(0),
            ProgressDrawTarget::term_like_with_hz(Box::new(RTerm), 8),
        );
        bar.set_style(
            ProgressStyle::with_template(
                "{spinner} [{elapsed_precise}] [{bar:30}] {pos}/{len} ({eta})",
            )
            .unwrap()
            .progress_chars("=>-"),
        );
        let (tx, rx) = unbounded();
        (IndicatifProgress { tx }, Self { bar, rx })
    }

    pub fn apply(&self, msg: ProgressMsg) {
        match msg {
            ProgressMsg::SetTotal(n) => self.bar.set_length(n),
            ProgressMsg::Increment(n) => self.bar.inc(n),
        }
    }

    pub fn drain(&self) {
        while let Ok(msg) = self.rx.try_recv() {
            self.apply(msg);
        }
    }

    pub fn recv(&self) -> Result<ProgressMsg, RecvTimeoutError> {
        self.rx.recv_timeout(PUMP_TICK)
    }

    pub fn finish(&self) {
        self.drain();
        self.bar.finish_and_clear();
    }
}
