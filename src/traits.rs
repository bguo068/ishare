use std::cmp::Ordering;

pub trait TotalOrd: PartialOrd + PartialEq {
    fn total_cmp(&self, other: &Self) -> Ordering;
}

impl TotalOrd for f32 {
    fn total_cmp(&self, other: &Self) -> Ordering {
        self.total_cmp(other)
    }
}

impl TotalOrd for f64 {
    fn total_cmp(&self, other: &Self) -> Ordering {
        self.total_cmp(other)
    }
}

impl TotalOrd for u64 {
    fn total_cmp(&self, other: &Self) -> Ordering {
        self.cmp(other)
    }
}

impl TotalOrd for u32 {
    fn total_cmp(&self, other: &Self) -> Ordering {
        self.cmp(other)
    }
}

impl TotalOrd for i32 {
    fn total_cmp(&self, other: &Self) -> Ordering {
        self.cmp(other)
    }
}
impl TotalOrd for u8 {
    fn total_cmp(&self, other: &Self) -> Ordering {
        self.cmp(other)
    }
}

impl TotalOrd for usize {
    fn total_cmp(&self, other: &Self) -> Ordering {
        self.cmp(other)
    }
}

impl<T1, T2> TotalOrd for (T1, T2)
where
    T1: TotalOrd,
    T2: TotalOrd,
{
    fn total_cmp(&self, other: &Self) -> Ordering {
        match self.0.total_cmp(&other.0) {
            Ordering::Equal => self.1.total_cmp(&other.1),
            x => x,
        }
    }
}

impl<T1, T2, T3> TotalOrd for (T1, T2, T3)
where
    T1: TotalOrd,
    T2: TotalOrd,
    T3: TotalOrd,
{
    fn total_cmp(&self, other: &Self) -> Ordering {
        match self.0.total_cmp(&other.0) {
            Ordering::Equal => match self.1.total_cmp(&other.1) {
                Ordering::Equal => self.2.total_cmp(&other.2),
                x => x,
            },
            x => x,
        }
    }
}

impl<T1, T2, T3, T4> TotalOrd for (T1, T2, T3, T4)
where
    T1: TotalOrd,
    T2: TotalOrd,
    T3: TotalOrd,
    T4: TotalOrd,
{
    fn total_cmp(&self, other: &Self) -> Ordering {
        match self.0.total_cmp(&other.0) {
            Ordering::Equal => match self.1.total_cmp(&other.1) {
                Ordering::Equal => match self.2.total_cmp(&other.2) {
                    Ordering::Equal => self.3.total_cmp(&other.3),
                    x => x,
                },
                x => x,
            },
            x => x,
        }
    }
}
impl<T> TotalOrd for &T
where
    T: TotalOrd,
{
    fn total_cmp(&self, other: &Self) -> Ordering {
        (*self).total_cmp(*other)
    }
}
