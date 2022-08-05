use std::{
    fmt::Display,
    ops::{Add, Neg, Sub},
};

use num_traits::Float;
use super::{ vector3::Vector3};
use crate::traits::{TrigConsts, Positional};

#[cfg(serde)]
use serde::{Serialize, Deserialize};

/*********************
 * STRUCT DEFINITION *
 *********************/

#[cfg_attr(serde, derive(Serialize, Deserialize))]
#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd)]
/// A point in 3D space
pub struct Cylindrical<T: num_traits::Float> {
    /// Angle from the positive `x` direction
    pub azimuth: T,
    /// Distance from the origin along the `xy` plane
    pub radius: T,
    /// Distance along the `z` axis
    pub height: T,
}

/***************************
 * CRATE TRAIT DEFINITIONS *
 ***************************/

macro_rules! impl_3d {
    ($var: tt) => {
        impl super::ThreeDimensionalConsts<$var> for Cylindrical<$var> {
            const ORIGIN: Self = Cylindrical {
                azimuth: 0.0,
                radius: 0.0,
                height: 0.0,
            };

            const UP: Self = Cylindrical {
                azimuth: 0.0,
                radius: 0.0,
                height: 1.0,
            };

            const DOWN: Self = Cylindrical {
                azimuth: 0.0,
                radius: 0.0,
                height: -1.0,
            };

            const FORWARD: Self = Cylindrical {
                azimuth: $var::FRAC_PI_2,
                radius: 1.0,
                height: 0.0,
            };

            const BACK: Self = Cylindrical {
                azimuth: $var::FRAC_3PI_2,
                radius: 1.0,
                height: 0.0,
            };

            const LEFT: Self = Cylindrical {
                azimuth: $var::PI,
                radius: 1.0,
                height: 0.0,
            };

            const RIGHT: Self = Cylindrical {
                azimuth: 0.0,
                radius: 1.0,
                height: 0.0,
            };
        }

        
    };
    ($($var : ident),+) => {
        $(impl_3d!($var);)+
    }
}

impl_3d!(f32, f64);

impl<T: Float> crate::traits::Magnitude<T> for Cylindrical<T> {
    fn magnitude(&self) -> T {
        self.quick_magnitude().sqrt()
    }

    fn quick_magnitude(&self) -> T {
        self.radius * self.radius + self.height * self.height
    }
}

impl<T: Float> crate::traits::Dot<T> for Cylindrical<T> {
    fn dot(&self, other: &Self) -> T {
        Into::<Vector3<T>>::into(self).dot(&other.into())
    }
}

impl<T: Float + TrigConsts> crate::traits::Cross3D for Cylindrical<T>{
    fn cross(&self, other: &Self) -> Self {
        Into::<Vector3<T>>::into(self).cross(&other.into()).into()
    }
}

/************************************
 * ARITHMETIC TRAIT IMPLEMENTATIONS *
 ************************************/

impl<T: Float + TrigConsts> Positional<T> for Cylindrical<T> {}


impl<T: Float+ TrigConsts> Add for Cylindrical<T> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        (Into::<Vector3<T>>::into(self) + Into::<Vector3<T>>::into(rhs)).into()
    }
}

impl<T: Float+ TrigConsts> Sub for Cylindrical<T> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        (Into::<Vector3<T>>::into(self) - Into::<Vector3<T>>::into(rhs)).into()
    }
}

impl<T: Float> std::ops::Div<T> for Cylindrical<T> {
    type Output = Self;

    fn div(self, rhs: T) -> Self::Output {
        Self {
            radius: self.radius / rhs,
            height: self.height / rhs,
            azimuth: self.azimuth,
        }
    }
}

impl<T: Float + TrigConsts> Neg for Cylindrical<T> {
    type Output = Cylindrical<T>;

    fn neg(self) -> Self::Output {
        Cylindrical {
            radius: self.radius,
            azimuth: (self.azimuth + T::PI) % T::TAU,
            height: -self.height,
        }
    }
}

/************************
 * FROM IMPLEMENTATIONS *
 ************************/

impl<T: Float> From<(T, T, T)> for Cylindrical<T> {
    fn from(tuple: (T, T, T)) -> Self {
        Cylindrical {
            radius: tuple.0,
            azimuth: tuple.1,
            height: tuple.2,
        }
    }
}

impl<T: Float> Into<(T, T, T)> for Cylindrical<T> {
    fn into(self) -> (T, T, T) {
        (self.radius, self.azimuth, self.height)
    }
}

impl<T: Float + TrigConsts> From<Vector3<T>> for Cylindrical<T> {
    fn from(cart: Vector3<T>) -> Self {
        let radius = (cart.x * cart.x + cart.y * cart.y).sqrt();
        let azimuthal_angle = if radius == T::ZERO {
            T::ZERO
        } else {
            if cart.x >= T::ZERO {
                (cart.y/radius).asin()
            } else {
                -(cart.y/radius).asin() + T::PI
            }
        };

        Cylindrical { azimuth: azimuthal_angle, radius, height: cart.z }
    }
}


impl<T: Float + TrigConsts> From<&Vector3<T>> for Cylindrical<T> {
    fn from(cart: &Vector3<T>) -> Self {
        let radius = (cart.x * cart.x + cart.y * cart.y).sqrt();
        let azimuthal_angle = if radius == T::ZERO {
            T::ZERO
        } else {
            if cart.x >= T::ZERO {
                (cart.y/radius).asin()
            } else {
                -(cart.y/radius).asin() + T::PI
            }
        };

        Cylindrical { azimuth: azimuthal_angle, radius, height: cart.z }
    }
}

/**************************
 * DISPLAY IMPLEMENTATION *
 **************************/

impl<T: Float + Display> Display for Cylindrical<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "({}, {}, {})", self.radius, self.azimuth, self.height)
    }
}



#[cfg(test)]
mod tests {
    use crate::traits::TrigConsts;
    use crate::traits::Positional;
    use crate::three_dimensional::ThreeDimensionalConsts;

    use super::Vector3;
    use super::Cylindrical;

    use std::f32::EPSILON as EPSILON;
    use assert_float_eq::*;

    #[test]
    pub fn convert_to_cartesian() {
        let cylindricals = [Cylindrical::<f32>::UP, Cylindrical::<f32>::DOWN, Cylindrical::<f32>::BACK, Cylindrical::<f32>::FORWARD, Cylindrical::<f32>::LEFT, Cylindrical::<f32>::RIGHT];
        let cartesians = [Vector3::<f32>::UP, Vector3::<f32>::DOWN, Vector3::<f32>::BACK, Vector3::<f32>::FORWARD, Vector3::<f32>::LEFT, Vector3::<f32>::RIGHT];

        for (s, c) in cylindricals.into_iter().zip(cartesians.into_iter()) {
            println!("{:?}", s);
            let s_star : Vector3<f32> = (&s).into();
            println!("{:?}", s_star);
            

            assert_float_absolute_eq!(c.x, s_star.x, EPSILON * s.azimuth.abs().log(2.0).max(1.0));
            assert_float_absolute_eq!(c.y, s_star.y, EPSILON * s.azimuth.abs().log(2.0).max(1.0));
            assert_float_absolute_eq!(c.z, s_star.z, EPSILON * s.height.abs().log(2.0).max(1.0));
        } 
    }

    #[test]
    pub fn is_positional() {
        let up = Cylindrical::<f32>::UP;

        for point in [Cylindrical::<f32>::BACK, Cylindrical::<f32>::FORWARD, Cylindrical::<f32>::LEFT, Cylindrical::<f32>::RIGHT] {
            println!("Angle between\n{:?} and\n{:?} is\n{}", &point, &up, up.angle_between(&point));
            assert_float_relative_eq!(f32::FRAC_PI_2, up.angle_between(&point), EPSILON);
        }

        assert_float_relative_eq!(f32::PI, up.angle_between(&Cylindrical::<f32>::DOWN),  EPSILON);
    }
}
