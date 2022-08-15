use std::ops::{Add, Neg, Sub};

use num_traits::Float;
use crate::traits::Positional;

mod vector3;
mod cylindrical;
mod spherical;

pub use vector3::*;
pub use cylindrical::*;
pub use spherical::*;

pub trait FullThreeDimensional<U> {}

impl<T, U: Float> FullThreeDimensional<U> for T where
    T: Positional<U>
        + ThreeDimensionalConsts<U>
        + Add
        + Sub
        + Neg
        + Sized
{
}

pub trait ThreeDimensionalConsts<T: Float> {
    const ORIGIN: Self;
    const UP: Self;
    const DOWN: Self;
    const FORWARD: Self;
    const BACK: Self;
    const LEFT: Self;
    const RIGHT: Self;
}
