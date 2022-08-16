use std::ops::{Add, Neg, Sub};

use num_traits::Float;
use crate::traits::Positional;

mod vector3;
mod cylindrical;
mod spherical;

pub use vector3::*;
pub use cylindrical::*;
pub use spherical::*;

/// Auto-trait for structs that implement consts and operators.
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

/// Trait holding constants for unit vectors and the origin
pub trait ThreeDimensionalConsts<T: Float> {

    /// Center of the coordinate space
    const ORIGIN: Self;
    /// Unit vector pointing in the positive z direction
    const UP: Self;
    /// Unit vector pointing in the negative z direction
    const DOWN: Self;
    /// Unit vector pointing in the positive x direction
    const LEFT: Self;
    /// Unit vector pointing in the negative x direction
    const RIGHT: Self;
    /// Unit vector pointing in the positive y direction
    const FORWARD: Self;
    /// Unit vector pointing in the negative y direction
    const BACK: Self;
}
