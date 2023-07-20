pub trait M2ARole: sealed::Sealed {}

pub struct OTSender;
impl M2ARole for OTSender {}

pub struct OTReceiver;
impl M2ARole for OTReceiver {}

mod sealed {
    pub trait Sealed {}
    impl Sealed for super::OTSender {}
    impl Sealed for super::OTReceiver {}
}
