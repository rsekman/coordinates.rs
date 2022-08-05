use std::ops::{Add, Neg, Sub};

use num_traits::Float;
use crate::traits::Positional;

pub mod vector3;
pub mod cylindrical;
pub mod spherical;

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
