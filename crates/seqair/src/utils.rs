pub trait TraceErr<T> {
    /// Like [`Option::and`], but traces a message when the value is `None`.
    fn trace_err(self, message: &'static str) -> Option<T>;
}

impl<T> TraceErr<T> for Option<T> {
    #[track_caller]
    #[inline(always)]
    fn trace_err(self, message: &'static str) -> Option<T> {
        if self.is_none() {
            tracing::trace!("{message}");
        }
        self
    }
}

pub trait TraceOk<T> {
    /// Like [`Result::ok`], but traces the error before discarding it.
    fn trace_ok(self, message: &'static str) -> Option<T>;
}

impl<T, E: std::fmt::Display> TraceOk<T> for Result<T, E> {
    #[track_caller]
    #[inline(always)]
    fn trace_ok(self, message: &'static str) -> Option<T> {
        match self {
            Ok(v) => Some(v),
            Err(error) => {
                tracing::trace!(%error, "{message}");
                None
            }
        }
    }
}
