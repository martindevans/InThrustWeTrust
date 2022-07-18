/// Calculate the area ratio (area / throat_area)
///
/// # Arguments
/// - `mach`: mach number at exit
/// - `gamma`: ratio of specific heats of the gas
/// 
/// # Reference
/// - https://www.grc.nasa.gov/www/k-12/rocket/astar.html
/// - https://en.wikipedia.org/wiki/File:Isentropic_Flow_Relations_Table.PNG
#[no_mangle]
pub extern fn nozzle_area_ratio(mach: f64, gamma: f64) -> f64
{
    let gam = gamma;
    let M = mach;

    // A/A* = {[(gam+1)/2]^-[(gam+1)/(gam-1)/2]} / Me * [1 + Me^2 * (gam-1)/2]^[(gam+1)/(gam-1)/2] 

    // {[(gam+1)/2]^-[(gam+1)/(gam-1)/2]} / Me
    let left = f64::powf(
        (gam + 1.0) / 2.0,
        -((gam + 1.0) / (gam - 1.0) / 2.0)
    ) / M;

    // [1 + Me^2 * (gam-1)/2]^[(gam+1)/(gam-1)/2] 
    let right = f64::powf(
        (1.0 + f64::powf(M, 2.0) * (gam - 1.0) / 2.0),
        ((gam + 1.0)/(gam - 1.0) / 2.0)
    );

    return left * right;
}

#[no_mangle]
pub extern fn nozzle_exit_temperature(temp_throat: f64, gamma: f64, mach_exit: f64) -> f64
{
    let Tt = temp_throat;
    let y = gamma;
    let M = mach_exit;

    let inner = 1.0 + (y - 1.0) / 2.0 * M * M;

    return Tt * f64::powf(inner, -1.0);
}

/// Calculate the mass flow rate of the exhaust at the throat assuming choked flow
/// 
/// # Arguments
/// - `throat_area`: area of the nozzle at the throat
/// - `query_area`: area of the nozzle at the query location
/// - `gamma`: ratio of specific heats of the gas
/// - `r`: gas constant for gas
/// 
/// # Reference
/// - https://www.grc.nasa.gov/www/k-12/rocket/wcora.html
#[no_mangle]
pub extern fn nozzle_mass_flow(throat_area: f64, pressure: f64, temperature: f64, gamma: f64, r: f64) -> f64
{
    let A = throat_area;
    let Pt = pressure;
    let Tt = temperature;
    let y = gamma;
    let R = r;

    let mdot = (Pt / f64::sqrt(R * Tt))
             * A
             * f64::sqrt(y)
             * f64::powf(
                (1.0 + (y - 1.0) / 2.0),
                (2.0 - 2.0 * y)
             );

    return mdot;
}

/// Calculate the total thrust for a given mass flow and exit velocity
/// 
/// # Arguments
/// - `mass_flow`: mass flow (kg/s) through the nozzle
/// - `exit_velocity`: velocity (m/s) of the exhaust
#[no_mangle]
pub extern fn nozzle_thrust(mass_flow: f64, exit_velocity: f64) -> f64
{
    return mass_flow * exit_velocity;
}

/// Calculate the exit velocity of the exhaust from the nozzle
/// 
/// # Arguments
/// - `inlet_temperature`: Temperature (K) of the exhaust at the inlet (throat)
/// - `inlet_pressure`: Pressure (Pa) at the throat
/// - `exit_pressure`: Pressure (Pa) at the end of the nozzle
/// - `molecular_weight`: Molecular weight (kg/kmol) of gas
/// - `gamma`: Cp/Cv, isentropic expansion factor
#[no_mangle]
pub extern fn nozzle_exit_velocity(throat_temperature: f64, inlet_pressure: f64, exit_pressure: f64, molecular_weight: f64, gamma: f64) -> f64
{
    // https://en.wikipedia.org/wiki/Rocket_engine_nozzle#de_Laval_nozzle_in_1_dimension

    let T = throat_temperature;
    let R = 8314.5;
    let M = molecular_weight;
    let y = gamma;
    let pe = exit_pressure;
    let p = inlet_pressure;

    let pow_term = f64::powf(
        pe / p,
        (y - 1.0) / y
    );

    return f64::sqrt(
          ((T * R) / M)
        * ((2.0 * y) / (y - 1.0))
        * (1.0 - pow_term)
    );
}

/// Calculate exit pressure of exhaust from the nozzle
/// 
/// # Arguments
/// - `throat_pressure`: Pressure (Pa) at the throat
/// - `gamma`: Cp/Cv, isentropic expansion factor
/// - `mach_exit`: mach of exhaust at nozzle exit
#[no_mangle]
pub extern fn nozzle_exit_pressure(throat_pressure: f64, gamma: f64, mach_exit: f64) -> f64
{
    let Pt = throat_pressure;
    let y = gamma;
    let M = mach_exit;

    let inner = 1.0 + (y - 1.0) / 2.0 * M * M;
    let exponent = -y / (y - 1.0);

    return Pt * f64::powf(inner, exponent);
}

#[no_mangle]
pub extern fn nozzle_exit_mach(area_ratio: f64, gamma: f64) -> f64
{
    let mut upper_limit = 1.0;
    while nozzle_area_ratio(upper_limit, gamma) < area_ratio {
        upper_limit += 1.0;
        if upper_limit > 25.0 {
            return -1.0;
        }
    }

    //let m1 = super::nozzle_area_ratio(10.0, 1.22);
    //panic!("{} {}", m0, m1);

    let mut conv = roots::SimpleConvergency { eps:0.01, max_iter:1000 };
    let root = roots::find_root_brent(
        1.0,
        upper_limit + 1.0,
        |x: f64| { nozzle_area_ratio(x, gamma) - area_ratio },
        &mut conv
    );

    return root.unwrap_or(-1.0);
}

#[cfg(test)]
mod tests
{
    #[test]
    fn nozzle_choked_flow_rate_raptor()
    {
        // Approximate numbers for the Raptor engine
        assert_eq!(
            f64::round(super::nozzle_mass_flow(0.2216, 3000000.0, 3587.0, 1.4, 287.06)),
            670.0
        );
    }

    #[test]
    fn area_ratio_from_mach_number_solver()
    {
        // Numbers taken from solver on this page: https://www.grc.nasa.gov/www/k-12/airplane/rktthsum.html
        assert_eq!(
            f64::round(super::nozzle_area_ratio(2.940, 1.4)),
            4.0
        );
    }

    #[test]
    fn nozzle_exit_velocity_raptor()
    {
        // Numbers taken from example on this page: https://en.wikipedia.org/wiki/Rocket_engine_nozzle#de_Laval_nozzle_in_1_dimension
        assert_eq!(
            f64::round(super::nozzle_exit_velocity(
                3500.0,
                7000000.0,
                100000.0,
                22.0,
                1.22
            )),
            2802.0
        );
    }

    #[test]
    fn nozzle_area_ratio_calculator()
    {
        // Numbers taken from calculator on this page: https://www.grc.nasa.gov/www/k-12/airplane/rktthsum.html
        assert_eq!(
            f64::round(super::nozzle_area_ratio(
                3.0,
                1.22
            ) * 1000.0),
            6339.0
        );
    }

    #[test]
    fn nozzle_exit_temperature_calculator()
    {
        // Numbers taken from calculator on this page: https://www.grc.nasa.gov/www/k-12/airplane/rktthsum.html
        assert_eq!(
            f64::round(super::nozzle_exit_temperature(
                3500.0,
                1.22,
                3.0
            )),
            1759.0
        );
    }

    #[test]
    fn nozzle_exit_pressure_calculator()
    {
        // Numbers taken from calculator on this page: https://www.grc.nasa.gov/www/k-12/airplane/rktthsum.html
        assert_eq!(
            f64::round(super::nozzle_exit_pressure(
                7000000.0,
                1.22,
                3.0
            )),
            154107.0
        );
    }

    #[test]
    fn root_finder_mach()
    {
        let mach = super::nozzle_exit_mach(6.339, 1.22);
        assert_eq!(f64::round(mach * 1000.0), 3000.0);
    }
}