use std::{
    fmt::Display,
    ops::{Add, Neg, Sub},
};

use num_traits::Float;

use super::{vector3::Vector3, cylindrical::Cylindrical};
use crate::traits::{Magnitude, Positional, TrigConsts};

#[cfg(serde)]
use serde::{Deserialize, Serialize};

/*********************
 * STRUCT DEFINITION *
 *********************/

#[cfg_attr(serde, derive(Serialize, Deserialize))]
#[derive(Debug, Copy,  Clone, PartialEq, Eq, PartialOrd)]
/// A point in 3D space using spherical coordinates as defined by [ISO 80000-2:2019](https://en.wikipedia.org/wiki/Spherical_coordinate_system#Definition).
pub struct Spherical<T: Float> {
    /// Distance from the origin
    pub radius: T,
    /// Angle from the positive `z` direction
    pub polar_angle: T,
    /// angle from the positive `x` direction along the `xy` plane
    pub azimuthal_angle: T,
}

impl<T: Float + TrigConsts> Spherical<T> {
    /// Maps parameters into appropriate domains.
    ///
    /// # Mapping
    ///
    /// - `Polar angle ∈ [0,π]`
    ///   - If `polar angle < 0`
    ///     - `Polar angle := -Polar angle`
    ///     - `Azimuthal angle := azimuthal angle + π`
    ///   - If `Polar angle > π`
    ///     - `Polar angle := polar angle % τ`
    ///     - If polar angle is still greater than π
    ///         - `Azimuthal angle := Azimuthal angle + π`
    ///         - `Polar angle := τ - Polar angle`
    /// - `azimuthal angle ∈ [0,τ]`
    ///   - If `azimuthal angle > τ`
    ///     - `azimuthal angle := azimuthal angle % τ`
    ///   - If `azimuthal angle < 0`
    ///     - `azimuthal angle := τ + azimuthal angle % τ`
    ///     - note: this is because rust by default will return the
    /// negative of the complement of the result given when taking
    /// the modulo with a positive numerator. i.e. `-2 % 5 = -3` whereas `2 % 5 = 2`
    /// - `radius ∈ [0,∞)`
    ///   - If `radius < 0`
    ///     - `radius := |radius|`
    ///     - `azimuthal angle := azimuthal angle + π % τ`
    ///     - `Polar angle := π - polar angle`
    pub fn new(radius: T, polar_angle: T, azimuthal_angle: T) -> Spherical<T> {
        // `Checked polar angle` ∈ [0,tau) when `polar angle` >= 0
        // `Checked polar angle` ∈ (0,tau] when `polar angle` <= -0
        let mut checked_polar_angle = polar_angle % T::TAU
            + if polar_angle.is_sign_negative() {
                // If polar angle is negative, checked polar angle ∈ (-tau, 0]
                // Add tau to move to the range (0,tau]
                T::TAU
            } else {
                T::ZERO
            };

        // `Checked azimuthal angle` ∈ [0,tau) when `azimuthal angle` >= 0
        // `Checked azimuthal angle` ∈ (0,tau] when `azimuthal angle` <= -0
        let mut checked_azimuthal_angle = (azimuthal_angle
            + if checked_polar_angle > T::PI {
                // `checked polar angle` is now ∈ [0,pi]
                checked_polar_angle = checked_polar_angle - T::PI;
                T::PI
            } else {
                T::ZERO
            })
            % T::TAU;

        let checked_radius = if radius.is_sign_negative() {
            checked_polar_angle = T::PI - checked_polar_angle;
            checked_azimuthal_angle = (checked_azimuthal_angle + T::PI) % T::TAU;
            -radius
        } else {
            radius
        };

        Spherical {
            polar_angle: checked_polar_angle,
            azimuthal_angle: checked_azimuthal_angle,
            radius: checked_radius,
        }
    }

    pub fn get_elevation(&self) -> T {
        T::FRAC_PI_2 - self.polar_angle
    }

    pub fn set_azimuthal_angle(&mut self, azimuthal_angle: T) {
        self.azimuthal_angle = azimuthal_angle % T::TAU;
    }

    pub fn set_polar_angle(&mut self, polar_angle: T) {
        // `Checked polar angle` ∈ [0,tau) when `polar angle` >= 0
        // `Checked polar angle` ∈ (0,tau] when `polar angle` <= -0
        let mut checked_polar_angle = polar_angle % T::TAU
            + if polar_angle.is_sign_negative() {
                // If polar angle is negative, checked polar angle ∈ (-tau, 0]
                // Add tau to move to the range (0,tau]
                T::TAU
            } else {
                T::ZERO
            };

        if checked_polar_angle > T::PI {
            // Checked polar angle is out of range, so azimuthal angle
            // will need to be adjusted as well
            self.set_azimuthal_angle(self.azimuthal_angle + T::PI);
            checked_polar_angle = T::TAU - checked_polar_angle;
        }

        self.polar_angle = checked_polar_angle;
    }
}

impl<T: Float + TrigConsts> super::ThreeDimensionalConsts<T> for Spherical<T> {
    const ORIGIN: Self = Spherical {
        azimuthal_angle: T::ZERO,
        radius: T::ZERO,
        polar_angle: T::ZERO,
    };

    const UP: Self = Spherical {
        azimuthal_angle: T::ZERO,
        radius: T::ONE,
        polar_angle: T::ZERO,
    };

    const DOWN: Self = Spherical {
        azimuthal_angle: T::ZERO,
        radius: T::ONE,
        polar_angle: T::PI,
    };

    const FORWARD: Self = Spherical {
        azimuthal_angle: T::FRAC_PI_2,
        radius: T::ONE,
        polar_angle: T::FRAC_PI_2,
    };

    const BACK: Self = Spherical {
        azimuthal_angle: T::FRAC_3PI_2,
        radius: T::ONE,
        polar_angle: T::FRAC_PI_2,
    };

    const LEFT: Self = Spherical {
        azimuthal_angle: T::PI,
        radius: T::ONE,
        polar_angle: T::FRAC_PI_2,
    };

    const RIGHT: Self = Spherical {
        azimuthal_angle: T::ZERO,
        radius: T::ONE,
        polar_angle: T::FRAC_PI_2,
    };
}

impl<T: Float> crate::traits::Magnitude<T> for Spherical<T> {
    #[inline]
    fn magnitude(&self) -> T {
        self.quick_magnitude()
    }

    #[inline]
    fn quick_magnitude(&self) -> T {
        self.radius
    }
}

impl<T: Float> crate::traits::Dot<T> for Spherical<T> {
    fn dot(&self, other: &Self) -> T {
        Into::<Vector3<T>>::into(self).dot(&other.into())
    }
}

impl<T: Float> crate::traits::Cross3D for Spherical<T> {
    fn cross(&self, other: &Self) -> Self {
        Into::<Vector3<T>>::into(self).cross(&other.into()).into()
    }
}

impl<T: Float + TrigConsts> Positional<T> for Spherical<T> {
    /// Using the spherical law of cosines
    fn angle_between(&self, other: &Self) -> T {
        let (lat_sin, lat_cos) = self.polar_angle.sin_cos();
        let (lat_sin_b, lat_cos_b) = other.polar_angle.sin_cos();

        (lat_cos * lat_cos_b
            + lat_sin * lat_sin_b * (self.azimuthal_angle - other.azimuthal_angle).cos())
        .acos()
    }
}

impl<T: Float + TrigConsts> Neg for Spherical<T> {
    type Output = Spherical<T>;

    fn neg(self) -> Self::Output {
        Spherical {
            polar_angle: T::PI - self.polar_angle,
            azimuthal_angle: (self.azimuthal_angle + T::PI) % T::TAU,
            radius: self.radius,
        }
    }
}

impl<T: Float> Add for Spherical<T> {
    type Output = Spherical<T>;

    fn add(self, rhs: Self) -> Self::Output {
        (Into::<Vector3<T>>::into(self) + Into::<Vector3<T>>::into(rhs)).into()
    }
}

impl<T: Float> Sub for Spherical<T> {
    type Output = Spherical<T>;

    fn sub(self, rhs: Self) -> Self::Output {
        (Into::<Vector3<T>>::into(self) - Into::<Vector3<T>>::into(rhs)).into()
    }
}

impl<T: Float> std::ops::Div<T> for Spherical<T> {
    type Output = Self;

    fn div(self, rhs: T) -> Self::Output {
        Self {
            radius: self.radius / rhs,
            azimuthal_angle: self.azimuthal_angle,
            polar_angle: self.polar_angle,
        }
    }
}

impl<T: Float> From<Vector3<T>> for Spherical<T> {
    fn from(cart: Vector3<T>) -> Self {
        From::from(&cart)
    }
}

impl<T: Float> From<&Vector3<T>> for Spherical<T> {
    fn from(cart: &Vector3<T>) -> Self {
        let radius = cart.magnitude();
        Spherical {
            polar_angle: (cart.z / radius).acos(),
            azimuthal_angle: cart.y.atan2(cart.x),
            radius,
        }
    }
}

impl<T: Float> From<Cylindrical<T>> for Spherical<T> {
    fn from(cyl: Cylindrical<T>) -> Self {
        From::<&Cylindrical<T>>::from(&cyl)
    }
}

impl<T: Float> From<&Cylindrical<T>> for Spherical<T> {
    fn from(cyl: &Cylindrical<T>) -> Self {
        Spherical {
            radius: cyl.magnitude(),
            polar_angle: cyl.radius.atan2(cyl.height),
            azimuthal_angle: cyl.azimuth,
        }
    }
}

impl<T: Float> From<(T, T, T)> for Spherical<T> {
    fn from(tuple: (T, T, T)) -> Self {
        Spherical {
            radius: tuple.0,
            polar_angle: tuple.1,
            azimuthal_angle: tuple.2,
        }
    }
}

impl<T: Display + Float> Display for Spherical<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "({},{},{})",
            self.radius, self.polar_angle, self.azimuthal_angle
        )
    }
}


#[cfg(test)]
mod tests {
    use crate::three_dimensional::ThreeDimensionalConsts;
    use crate::traits::TrigConsts;
    use crate::traits::Positional;

    use super::Spherical;
    use super::Vector3;

    use std::f32::EPSILON as EPSILON;
    use assert_float_eq::*;

    #[test]
    pub fn convert_to_cartesian() {
        let sphericals = [Spherical::<f32>::UP, Spherical::<f32>::DOWN, Spherical::<f32>::BACK, Spherical::<f32>::FORWARD, Spherical::<f32>::LEFT, Spherical::<f32>::RIGHT];
        let cartesians = [Vector3::<f32>::UP, Vector3::<f32>::DOWN, Vector3::<f32>::BACK, Vector3::<f32>::FORWARD, Vector3::<f32>::LEFT, Vector3::<f32>::RIGHT];

        for (s, c) in sphericals.into_iter().zip(cartesians.into_iter()) {
            println!("{:?}", s);
            let s_star : Vector3<f32> = (&s).into();
            println!("{:?}", s_star);

            assert_float_absolute_eq!(c.x, s_star.x, EPSILON * s.azimuthal_angle.abs().log(2.0).max(1.0));
            assert_float_absolute_eq!(c.y, s_star.y, EPSILON * s.azimuthal_angle.abs().log(2.0).max(1.0));
            assert_float_absolute_eq!(c.z, s_star.z, EPSILON * s.polar_angle.abs().log(2.0).max(1.0));
        }  
    }

    #[test]
    pub fn is_positional() {
        let up = Spherical::<f32>::UP;

        for point in [Spherical::<f32>::BACK, Spherical::<f32>::FORWARD, Spherical::<f32>::LEFT, Spherical::<f32>::RIGHT] {
            assert_float_relative_eq!(f32::FRAC_PI_2, up.angle_between(&point), EPSILON);
        }

        assert_float_relative_eq!(f32::PI, up.angle_between(&Spherical::<f32>::DOWN),  EPSILON);
    }
}