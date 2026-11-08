//! Micro-library for converting from geodetic to UTM coordinates.
#![cfg_attr(feature = "no_std", no_std)]

#[cfg(feature = "no_std")]
extern crate core as std;

use std::f64::consts::PI;

#[cfg(feature = "no_std")]
extern crate num;

#[cfg(feature = "no_std")]
#[allow(unused_imports)]
// it's not clear why this generates an unused imports, b/c tests fail w/o it
use num::traits::float::Float;

pub struct Ellipsoid {
    a: f64,
    f: f64,
}

const WGS84: Ellipsoid = Ellipsoid {
    a: 6378137.0,
    f: 1.0 / 298.257223563,
};

const ZONE_LETTERS: &'static str = "CDEFGHJKLMNPQRSTUVWXX";

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

/// Converts a latitude and longitude in decimal degrees (using the latitude and longitude to compute the zone) to UTM coordinates using the WGS84 ellipsoid.
///
/// # Examples
///
/// ```
/// use utm::to_utm_wgs84_no_zone;
/// let (northing, easting, meridian_convergence) = to_utm_wgs84_no_zone(40.62, -123.45);
/// ```
pub fn to_utm_wgs84_no_zone(latitude: f64, longitude: f64) -> (f64, f64, f64) {
    to_utm_wgs84(
        latitude,
        longitude,
        lat_lon_to_zone_number(latitude, longitude),
    )
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
    // longitude from the central meridian, in [-PI, PI): zone 60 reaches over the date line
    let dlon = longitude - (long_origin * PI / 180.0);
    let dlon = dlon - 2.0 * PI * ((dlon + PI) / (2.0 * PI)).floor();
    let a = latitude.cos() * dlon;

    // meridian arc, series in e^2 up to the e^8 terms
    let e4 = e2 * e2;
    let term1 =
        1.0 - e2 / 4.0 - (3.0 * e4) / 64.0 - (5.0 * e4 * e2) / 256.0 - (175.0 * e4 * e4) / 16384.0;
    let term2 = (3.0 * e2) / 8.0
        + (3.0 * e4) / 32.0
        + (45.0 * e4 * e2) / 1024.0
        + (105.0 * e4 * e4) / 4096.0;
    let term3 = (15.0 * e4) / 256.0 + (45.0 * e4 * e2) / 1024.0 + (525.0 * e4 * e4) / 16384.0;
    let term4 = (35.0 * e4 * e2) / 3072.0 + (175.0 * e4 * e4) / 12288.0;
    let term5 = (315.0 * e4 * e4) / 131072.0;

    let m = ellipsoid.a
        * (term1 * latitude - term2 * (2.0 * latitude).sin() + term3 * (4.0 * latitude).sin()
            - term4 * (6.0 * latitude).sin()
            + term5 * (8.0 * latitude).sin());

    let x1 = ((1.0 - t + c) * a * a * a) / 6.0;
    let x2 = ((5.0 - 18.0 * t + t * t + 72.0 * c - 58.0 * ep2) * a * a * a * a * a) / 120.0;
    let x = 0.9996 * n * (a + x1 + x2);

    let y1 = (5.0 - t + 9.0 * c + 4.0 * c * c) * (a * a * a * a) / 24.0;
    let y2 = (61.0 - 58.0 * t + t * t + 600.0 * c - 330.0 * ep2) * (a * a * a * a * a * a) / 720.0;
    let y3 = (a * a) / 2.0 + y1 + y2;
    let y = 0.9996 * (m + n * latitude.tan() * y3);

    // Southern hemisphere coords add a factor. The equator is in the northern hemisphere, as for
    // `lat_to_zone_letter`, which gives 'N' for a latitude of 0.
    let northing_mod = if latitude >= 0.0 { 0.0 } else { 10000000.0 };

    let northing = y + northing_mod;
    let easting = x + 500000.0;

    // the convergence is computed from the northing measured from the equator, without the
    // false northing of the southern hemisphere
    let meridian_convergence = meridian_convergence(y, easting, WGS84);
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

const K0: f64 = 0.9996;

#[cfg(feature = "std")]
impl std::error::Error for WSG84ToLatLonError {
    fn description(&self) -> &str {
        match self {
            WSG84ToLatLonError::EastingOutOfRange => {
                "Easting out of range, must be between 100000 and 999999"
            }
            WSG84ToLatLonError::NorthingOutOfRange => {
                "Northing out of range, must be between 0 and 10000000"
            }
            WSG84ToLatLonError::ZoneNumOutOfRange => {
                "Zone num out of range, must be between 1 and 60"
            }
            WSG84ToLatLonError::ZoneLetterOutOfRange => {
                "Zone letter out of range, must be between C and X"
            }
        }
    }
}

#[derive(Debug, PartialEq)]
/// Error type for the wsg84_utm_to_lat_lon function.
pub enum WSG84ToLatLonError {
    EastingOutOfRange,
    NorthingOutOfRange,
    ZoneNumOutOfRange,
    ZoneLetterOutOfRange,
}

/// Converts a UTM coordinate to a latitude and longitude.
/// zone_num can be obtain by calling lat_lon_to_zone_number
/// zone_letter can be obtain by calling lat_to_zone_letter
///
/// # Examples
///
/// ```
/// use utm::wsg84_utm_to_lat_lon;
/// const DELTA: f64 = 3e-5;
/// fn is_close(a: f64, b: f64, epsilon: f64) -> bool {
///        (a - b).abs() < epsilon
/// }
/// // Capetown, South Africa,
/// let easting = 261878_f64;
/// let northing = 6243186_f64;
/// let zone_num = 34_u8;
/// let zone_letter = 'H';
/// let (lat, long) = wsg84_utm_to_lat_lon(easting, northing, zone_num, zone_letter).unwrap();
/// assert_eq!(is_close(lat, -33.92487, DELTA), true);
/// assert_eq!(is_close(long, 18.42406, DELTA), true);
/// ```
pub fn wsg84_utm_to_lat_lon(
    easting: f64,
    northing: f64,
    zone_num: u8,
    zone_letter: char,
) -> Result<(f64, f64), WSG84ToLatLonError> {
    if easting < 100000. || 1000000. <= easting {
        return Err(WSG84ToLatLonError::EastingOutOfRange);
    }
    if northing < 0. || northing > 10000000. {
        return Err(WSG84ToLatLonError::NorthingOutOfRange);
    }
    if zone_num < 1 || zone_num > 60 {
        return Err(WSG84ToLatLonError::ZoneNumOutOfRange);
    }
    if zone_letter < 'C' || zone_letter > 'X' {
        return Err(WSG84ToLatLonError::ZoneLetterOutOfRange);
    }

    let ellipsoid = WGS84;

    // first eccentricity squared, as in the forward conversion
    let e_sq = 2.0 * ellipsoid.f - ellipsoid.f * ellipsoid.f;
    let e2 = e_sq.powi(2);
    let e3 = e_sq.powi(3);
    let e_p2: f64 = e_sq / (1. - e_sq);

    let sqrt_e: f64 = (1. - e_sq).sqrt();
    let _e: f64 = (1. - sqrt_e) / (1. + sqrt_e);
    let _e2: f64 = _e.powi(2);
    let _e3: f64 = _e.powi(3);
    let _e4: f64 = _e.powi(4);
    let _e5: f64 = _e.powi(5);

    let m1 = 1. - e_sq / 4. - 3. * e2 / 64. - 5. * e3 / 256. - 175. * e2 * e2 / 16384.;

    let p2: f64 = 3. / 2. * _e - 27. / 32. * _e3 + 269. / 512. * _e5;
    let p3: f64 = 21. / 16. * _e2 - 55. / 32. * _e4;
    let p4: f64 = 151. / 96. * _e3 - 417. / 128. * _e5;
    let p5: f64 = 1097. / 512. * _e4;

    let x = easting - 500000_f64;
    let mut y = northing;

    let northern = zone_letter >= 'N';

    if !northern {
        y -= 1e7;
    }

    let m = y / K0;
    let mu = m / (ellipsoid.a * m1);

    let p_rad = mu
        + p2 * (2. * mu).sin()
        + p3 * (4. * mu).sin()
        + p4 * (6. * mu).sin()
        + p5 * (8. * mu).sin();

    let p_sin = p_rad.sin();
    let p_sin2 = p_sin.powi(2);

    let p_cos = p_rad.cos();

    let p_tan = p_rad.tan();
    let p_tan2 = p_tan.powi(2);
    let p_tan4 = p_tan.powi(4);

    let ep_sin = 1. - e_sq * p_sin2;
    let ep_sin_sqrt = ep_sin.sqrt();

    let n = ellipsoid.a / ep_sin_sqrt;
    let r = (1. - e_sq) / ep_sin;

    let c = e_p2 * p_cos * p_cos;
    let c2 = c * c;

    let d = x / (n * K0);
    let d2 = d.powi(2);
    let d3 = d.powi(3);
    let d4 = d.powi(4);
    let d5 = d.powi(5);
    let d6 = d.powi(6);

    let latitude = p_rad
        - (p_tan / r) * (d2 / 2. - d4 / 24. * (5. + 3. * p_tan2 + 10. * c - 4. * c2 - 9. * e_p2))
        + d6 / 720. * (61. + 90. * p_tan2 + 298. * c + 45. * p_tan4 - 252. * e_p2 - 3. * c2);

    let longitude = (d - d3 / 6. * (1. + 2. * p_tan2 + c)
        + d5 / 120. * (5. - 2. * c + 28. * p_tan2 - 3. * c2 + 8. * e_p2 + 24. * p_tan4))
        / p_cos;

    let longitude = longitude / PI * 180. + (f64::from(zone_num) - 1.) * 6. - 180. + 3.;
    Ok((
        latitude / PI * 180.,
        // zone 60 reaches over the date line: report 181.5 as -178.5
        longitude - 360. * ((longitude + 180.) / 360.).floor(),
    ))
}

/// Convert a latitude to the UTM zone letter.
///
/// # Examples
///
/// ```
/// use utm::lat_to_zone_letter;
/// assert_eq!(lat_to_zone_letter(-33.92487), Some('H'));
/// assert_eq!(lat_to_zone_letter(0.), Some('N'));
/// assert_eq!(lat_to_zone_letter(50.77535), Some('U'));
/// ```
pub fn lat_to_zone_letter(latitude: f64) -> Option<char> {
    if -80. <= latitude && latitude <= 84. {
        let (_, c) = ZONE_LETTERS
            .char_indices()
            .nth(((latitude + 80.) / 8.).floor() as usize)
            .unwrap();
        return Some(c);
    }
    return None;
}

/// Convert a latitude and longitude to the UTM zone number.
///
/// # Examples
///
/// ```
/// use utm::lat_lon_to_zone_number;
/// assert_eq!(lat_lon_to_zone_number(-33.92487, 18.42406), 34);
/// assert_eq!(lat_lon_to_zone_number(0., 0.), 31);
/// assert_eq!(lat_lon_to_zone_number(50.77535, 6.08389), 32);
/// // the date line is where zone 60 ends and zone 1 begins
/// assert_eq!(lat_lon_to_zone_number(0., 179.999), 60);
/// assert_eq!(lat_lon_to_zone_number(0., 180.), 1);
/// ```
pub fn lat_lon_to_zone_number(latitude: f64, longitude: f64) -> u8 {
    // 180 is the same meridian as -180, and 190 as -170: zone numbers are 1 to 60
    let longitude = longitude - 360. * ((longitude + 180.) / 360.).floor();

    if 56. <= latitude && latitude < 64. && 3. <= longitude && longitude < 12. {
        return 32;
    }

    if 72. <= latitude && latitude <= 84. && longitude >= 0. {
        if longitude < 9. {
            return 31;
        }
        if longitude < 21. {
            return 33;
        }
        if longitude < 33. {
            return 35;
        }
        if longitude < 42. {
            return 37;
        }
    }

    return (((longitude + 180.) / 6.).floor() + 1.) as u8;
}

#[cfg(test)]
mod tests {
    use super::*;

    const DELTA: f64 = 3e-5;

    #[test]
    fn reference() {
        let latitude = 60.9679875497;
        let longitude = -149.119325194;
        let (northing, easting, meridian_convergence) = to_utm_wgs84(latitude, longitude, 6);
        assert!((385273.02 - easting).abs() < 1e-2);
        assert!((6761077.20 - northing).abs() < 1e-2);
        assert!((0.0323 - meridian_convergence).abs() < 1e-4);
    }

    #[test]
    fn test_southern_hemisphere() {
        let latitude = -23.809749;
        let longitude = 148.334111;
        let (northing, easting, meridian_convergence) = to_utm_wgs84(latitude, longitude, 55);
        assert!((635900.94 - easting).abs() < 1e-2);
        assert!((7366197.82 - northing).abs() < 1e-2);
        // PROJ: convergence -0.53866478 degrees, which is +0.00940 here (the sign convention is the opposite of PROJ's)
        assert!((0.0094 - meridian_convergence).abs() < 1e-4);
    }

    #[test]
    fn test_to_lat_lon() {
        let (expected_lat, expected_lon) = (-41.28646, 174.77624);

        let wrong_easting = 50.;
        let wrong_northing = -1.;
        let wrong_zone_num = 61;
        let wrong_zone_letter = 'y';

        let easting = 313784.;
        let northing = 5427057.;
        let zone_num = 60;
        let zone_letter = 'G';

        let mut result = wsg84_utm_to_lat_lon(wrong_easting, northing, zone_num, zone_letter);
        assert!(result.is_err());
        assert_eq!(result.unwrap_err(), WSG84ToLatLonError::EastingOutOfRange);

        result = wsg84_utm_to_lat_lon(easting, wrong_northing, zone_num, zone_letter);
        assert!(result.is_err());
        assert_eq!(result.unwrap_err(), WSG84ToLatLonError::NorthingOutOfRange);

        result = wsg84_utm_to_lat_lon(easting, northing, wrong_zone_num, zone_letter);
        assert!(result.is_err());
        assert_eq!(result.unwrap_err(), WSG84ToLatLonError::ZoneNumOutOfRange);

        result = wsg84_utm_to_lat_lon(easting, northing, zone_num, wrong_zone_letter);
        assert!(result.is_err());
        assert_eq!(
            result.unwrap_err(),
            WSG84ToLatLonError::ZoneLetterOutOfRange
        );

        result = wsg84_utm_to_lat_lon(easting, northing, zone_num, zone_letter);
        assert!(result.is_ok());
        let (latitude, longitude) = result.unwrap();
        assert_eq!(is_close(latitude, expected_lat, DELTA), true);
        assert_eq!(is_close(longitude, expected_lon, DELTA), true);
    }

    #[test]
    fn test_to_wsg84_no_zone() {
        let latitude = 60.9679875497;
        let longitude = -149.119325194;
        let (northing, easting, meridian_convergence) = to_utm_wgs84(latitude, longitude, 6);

        let (northing_2, easting_2, meridian_convergence_2) =
            to_utm_wgs84_no_zone(latitude, longitude);
        assert_eq!(northing, northing_2);
        assert_eq!(easting, easting_2);
        assert_eq!(meridian_convergence, meridian_convergence_2);
    }

    // Reference values of the tests below come from PROJ 9.8.1:
    //   echo "<lon> <lat>" | cs2cs -f "%.6f" +proj=longlat +datum=WGS84 +to +proj=utm +zone=<zone> [+south] +datum=WGS84

    #[test]
    fn wgs84_ellipsoid() {
        // the defining parameters of WGS 84: a = 6378137 m, 1/f = 298.257223563 (not GRS 80's 298.257222101)
        assert_eq!(WGS84.a, 6378137.0);
        assert_eq!(WGS84.f, 1.0 / 298.257223563);
    }

    #[test]
    fn central_meridian_matches_proj() {
        // on the central meridian only the meridian arc counts, so this tests the ellipsoid and the arc series
        let reference = [
            (70.0, 7765873.135479, 35),
            (45.0, 4982950.400227, 35),
            (0.5, 55265.037143, 35),
            (-70.0, 2234126.864521, 35),
        ];
        for (latitude, northing, zone) in reference {
            let (n, e, _) = to_utm_wgs84(latitude, 27.0, zone);
            assert!((e - 500000.0).abs() < 1e-6, "easting {} at {}", e, latitude);
            assert!(
                (n - northing).abs() < 2e-4,
                "northing {} vs {} at {}",
                n,
                northing,
                latitude
            );
        }
    }

    #[test]
    fn to_lat_lon_matches_proj_at_the_edge_of_the_zone() {
        // 3 degrees from the central meridian of zone 35 (27 E): (lon, lat, easting, northing, zone letter)
        let reference = [
            (30.0, 60.0, 667294.821124, 6655205.483635, 'V'),
            (24.0, 70.0, 385526.285254, 7768690.108775, 'W'),
            (30.0, 45.0, 736446.026101, 4987329.504699, 'T'),
            (24.0, -60.0, 332705.178876, 3344794.516365, 'E'),
            (30.0, -30.0, 789409.653236, 6677424.095615, 'J'),
        ];
        for (lon, lat, easting, northing, letter) in reference {
            assert_eq!(lat_to_zone_letter(lat), Some(letter));
            let (la, lo) = wsg84_utm_to_lat_lon(easting, northing, 35, letter).unwrap();
            // 2e-7 degrees is about 2 cm
            assert!((la - lat).abs() < 2e-7, "latitude {} vs {}", la, lat);
            assert!((lo - lon).abs() < 2e-7, "longitude {} vs {}", lo, lon);
        }
    }

    #[test]
    fn round_trip() {
        // out and back again, in the southern hemisphere (the point of the test above) and in the northern
        for (latitude, longitude, zone) in [
            (-23.809749, 148.334111, 55),
            (60.9679875497, -149.119325194, 6),
        ] {
            let (northing, easting, _) = to_utm_wgs84(latitude, longitude, zone);
            let letter = lat_to_zone_letter(latitude).unwrap();
            let (la, lo) = wsg84_utm_to_lat_lon(easting, northing, zone, letter).unwrap();
            assert!(
                (la - latitude).abs() < 2e-8,
                "latitude {} vs {}",
                la,
                latitude
            );
            assert!(
                (lo - longitude).abs() < 2e-8,
                "longitude {} vs {}",
                lo,
                longitude
            );
        }
    }

    #[test]
    fn equator_is_in_the_northern_hemisphere() {
        // PROJ: 28 E, 0 N, zone 35 -> 611280.650891, 0
        for latitude in [0.0, -0.0] {
            let (northing, easting, _) = to_utm_wgs84(latitude, 28.0, 35);
            assert!(northing.abs() < 1e-6, "northing {}", northing);
            assert!(
                (easting - 611280.650891).abs() < 1e-3,
                "easting {}",
                easting
            );
        }
        // and the zone letter that goes with it takes the point back
        let letter = lat_to_zone_letter(0.0).unwrap();
        assert_eq!(letter, 'N');
        let (northing, easting, _) = to_utm_wgs84(0.0, 28.0, 35);
        let (latitude, longitude) = wsg84_utm_to_lat_lon(easting, northing, 35, letter).unwrap();
        assert!(latitude.abs() < 1e-9, "latitude {}", latitude);
        assert!((longitude - 28.0).abs() < 1e-9, "longitude {}", longitude);
    }

    #[test]
    fn zone_60_reaches_over_the_date_line() {
        // PROJ: -178.5 E, 23.8 S, zone 60 -> 958724.569875, 7360636.327472 (zone 60 is centred on 177 E)
        let (northing, easting, _) = to_utm_wgs84(-23.8, -178.5, 60);
        assert!(
            (easting - 958724.569875).abs() < 1e-3,
            "easting {}",
            easting
        );
        assert!(
            (northing - 7360636.327472).abs() < 1e-3,
            "northing {}",
            northing
        );
        // the same point, with the longitude as 181.5 E, and in zone 1 (centred on 177 W)
        let (northing, easting, _) = to_utm_wgs84(-23.8, 181.5, 60);
        assert!((easting - 958724.569875).abs() < 1e-3 && (northing - 7360636.327472).abs() < 1e-3);
        let (northing, easting, _) = to_utm_wgs84(-23.8, 178.5, 1);
        assert!((easting - 41275.430125).abs() < 1e-3, "easting {}", easting);
        assert!(
            (northing - 7360636.327472).abs() < 1e-3,
            "northing {}",
            northing
        );
        // back: a longitude in [-180, 180)
        let (latitude, longitude) =
            wsg84_utm_to_lat_lon(958724.569875, 7360636.327472, 60, 'K').unwrap();
        assert!(
            (latitude + 23.8).abs() < 1e-6 && (longitude + 178.5).abs() < 1e-6,
            "{} {}",
            latitude,
            longitude
        );
    }

    #[test]
    fn zone_numbers_are_1_to_60() {
        assert_eq!(lat_lon_to_zone_number(10., 180.), 1);
        assert_eq!(lat_lon_to_zone_number(10., -180.), 1);
        assert_eq!(lat_lon_to_zone_number(10., 179.999), 60);
        assert_eq!(lat_lon_to_zone_number(10., 190.), 2);
        assert_eq!(lat_lon_to_zone_number(10., -190.), 59); // -190 is 170 E
    }

    #[test]
    fn meridian_convergence_matches_proj() {
        // `proj -V +proj=utm +zone=<zone> [+south] +datum=WGS84` prints the convergence in degrees; this
        // function returns radians, with the opposite sign (positive west of the central meridian in the northern hemisphere)
        let reference = [
            (60.0, 30.0, 35, 2.59867269),
            (60.9679875497, -149.119325194, 6, -1.85322919),
            (-23.809749, 148.334111, 55, -0.53866478),
            (-60.0, 24.5, 35, 2.16540871),
        ];
        for (latitude, longitude, zone, proj_degrees) in reference {
            let (_, _, convergence) = to_utm_wgs84(latitude, longitude, zone);
            let expected = -f64::to_radians(proj_degrees);
            assert!(
                (convergence - expected).abs() < 1e-6,
                "{} vs {} at {}",
                convergence,
                expected,
                latitude
            );
        }
    }

    fn is_close(a: f64, b: f64, epsilon: f64) -> bool {
        (a - b).abs() < epsilon
    }
}
