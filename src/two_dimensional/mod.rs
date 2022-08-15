use std::ops::{Add, Sub, Neg};

use num_traits::Float;

use crate::traits::Positional;

mod polar;
mod vector2;

pub use polar::*;
pub use vector2::*;

pub trait FullTwoDimensional<U> {}

impl<T, U: Float> FullTwoDimensional<U> for T where
    T: Positional<U>
        + TwoDimensionalConsts<U>
        + Add
        + Sub
        + Neg
        + Sized
{
}



pub trait TwoDimensionalConsts<T: Float> {
    const ORIGIN: Self;
    const UP: Self;
    const DOWN: Self;
    const LEFT: Self;
    const RIGHT: Self;
}
