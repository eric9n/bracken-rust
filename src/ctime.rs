use std::time::{Duration, SystemTime};

pub fn timeval_subtract(start: SystemTime, end: SystemTime) -> Result<Duration, bool> {
    match end.duration_since(start) {
        Ok(duration) => Ok(duration),
        Err(_) => Err(true), // 返回true表示结果为负
    }
}
