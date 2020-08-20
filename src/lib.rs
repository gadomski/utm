//! Micro-library for converting from geodetic to UTM coordinates.
#![cfg_attr(feature = "no_std", no_std)]

#[cfg(feature = "no_std")]
extern crate core as std;

use std::f64::consts::PI;

#[cfg(feature = "no_std")]
extern crate num;


pub struct Ellipsoid {
    a: f64,
    f: f64,
}

const WGS84: Ellipsoid = Ellipsoid {
    a: 6378137.0,
    f: 1.0 / 298.257222101,
};

/// Converts a latitude and longitude in decimal degrees to UTM coordinates using the WGS84 ellipsoid.
///
/// # Examples
///
/// ```
/// use utm::to_utm_wgs84;
/// let (northing, easting, meridian_convergence) = to_utm_wgs84(40.62, -123.45, 10);
/// ```
pub fn to_utm_wgs84(latitude: f64, longitude: f64, zone: u8) -> (f64, f64, f64) {
    let latitude = latitude * PI / 180.0;
    let longitude = longitude * PI / 180.0;
    radians_to_utm_wgs84(latitude, longitude, zone)
}

/// Converts a latitude and longitude in radians to UTM coordinates using the WGS84 ellipsoid.
///
/// # Examples
///
/// ```
/// use std::f64::consts::PI;
/// use utm::radians_to_utm_wgs84;
/// let latitude = 40.62 * PI / 180.0;
/// let longitude = -123.45 * PI / 180.0;
/// let (northing, easting, meridian_convergence) = radians_to_utm_wgs84(latitude, longitude, 10);
/// ```
pub fn radians_to_utm_wgs84(latitude: f64, longitude: f64, zone: u8) -> (f64, f64, f64) {
    let ellipsoid = WGS84;
    let long_origin = zone as f64 * 6.0 - 183.0;
    let e2 = 2.0 * ellipsoid.f - ellipsoid.f * ellipsoid.f;
    let ep2 = e2 / (1.0 - e2);

    let n = ellipsoid.a / (1.0 - e2 * latitude.sin() * latitude.sin()).sqrt();
    let t = latitude.tan() * latitude.tan();
    let c = ep2 * latitude.cos() * latitude.cos();
    let a = latitude.cos() * (longitude - (long_origin * PI / 180.0));

    let term1 = 1.0 - e2 / 4.0 - (3.0 * e2 * e2) / 64.0 - (5.0 * e2 * e2 * e2) / 256.0;
    let term2 = (3.0 * e2) / 8.0 + (3.0 * e2 * e2) / 32.0 + (45.0 * e2 * e2 * e2) / 1024.0;
    let term3 = (15.0 * e2 * e2) / 256.0 + (45.0 * e2 * e2 * e2) / 1024.0;
    let term4 = (35.0 * e2 * e2 * e2) / 3072.0;

    let m = ellipsoid.a
        * (term1 * latitude - term2 * (2.0 * latitude).sin() + term3 * (4.0 * latitude).sin()
            - term4 * (6.0 * latitude).sin());

    let x1 = ((1.0 - t + c) * a * a * a) / 6.0;
    let x2 = ((5.0 - 18.0 * t + t * t + 72.0 * c - 58.0 * ep2) * a * a * a * a * a) / 120.0;
    let x = 0.9996 * n * (a + x1 + x2);

    let y1 = (5.0 - t + 9.0 * c + 4.0 * c * c) * (a * a * a * a) / 24.0;
    let y2 = (61.0 - 58.0 * t + t * t + 600.0 * c - 330.0 * ep2) * (a * a * a * a * a * a) / 720.0;
    let y3 = (a * a) / 2.0 + y1 + y2;
    let y = 0.9996 * (m + n * latitude.tan() * y3);

    let northing = y;
    let easting = x + 500000.0;

    let meridian_convergence = meridian_convergence(northing, easting, WGS84);
    (northing, easting, meridian_convergence)
}

fn meridian_convergence(northing: f64, easting: f64, ellipsoid: Ellipsoid) -> f64 {
    let e2: f64 = 2.0 * ellipsoid.f - ellipsoid.f * ellipsoid.f;
    let e1 = (1.0 - (1.0 - e2).sqrt()) / (1.0 + (1.0 - e2).sqrt());
    let mu_const =
        ellipsoid.a * (1.0 - e2 / 4.0 - 3.0 * e2 * e2 / 64.0 - 5.0 * e2 * e2 * e2 / 256.0);

    let np = northing / 0.9996;
    let mu = np / mu_const;
    let foot_lat = footprint_latitude(e1, mu);

    let ep = (easting - 500000.0) / 0.9996;
    let n = ellipsoid.a / (1.0 - e2 * foot_lat.sin() * foot_lat.sin()).sqrt();
    let m = (ellipsoid.a * (1.0 - e2)) / (1.0 - e2 * foot_lat.sin() * foot_lat.sin()).powf(1.5);

    let conv1 = -(ep / n) * foot_lat.tan();
    let h30 = (ep / n).powi(3);
    let k28 = n / m;
    let k29 = k28 * k28;
    let j29 = foot_lat.tan() * foot_lat.tan();
    let conv2 = (foot_lat.tan() * h30 / 3.0) * (-2.0 * k29 + 3.0 * k28 + j29);
    conv1 + conv2
}

fn footprint_latitude(e1: f64, mu: f64) -> f64 {
    let term1 = 3.0 * e1 / 2.0 - 27.0 * e1 * e1 * e1 / 32.0;
    let term2 = 21.0 * e1 * e1 / 16.0 - 55.0 * e1 * e1 * e1 * e1 / 32.0;
    let term3 = 151.0 * e1 * e1 * e1 / 96.0;
    let term4 = 1097.0 * e1 * e1 * e1 * e1 / 512.0;

    mu + term1 * (2.0 * mu).sin()
        + term2 * (4.0 * mu).sin()
        + term3 * (6.0 * mu).sin()
        + term4 * (8.0 * mu).sin()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reference() {
        let latitude = 60.9679875497;
        let longitude = -149.119325194;
        let (northing, easting, meridan_convgergence) = to_utm_wgs84(latitude, longitude, 6);
        assert!((385273.02 - easting).abs() < 1e-2);
        assert!((6761077.20 - northing).abs() < 1e-2);
        assert!((0.0323 - meridan_convgergence).abs() < 1e-4);
    }
}
