using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class OrbitPhysics : MonoBehaviour
{
    public const float GRAV_CONST = 6.6743e-11f; //physics
    private float standardGrav = -1;
    private OrbitGeometry geo;
 
    #region Math Helpers
    private static float Square(float val) => val * val; //dont trust the pow function to be effecient enough
    private static float Cube(float val) => val * Square(val);
    private static float Mod(float val, float r) => (val % r + r) % r; //get true Modulo function
    #endregion

    private OrbitPhysics() { }

    public OrbitPhysics(OrbitGeometry geo, float standardGrav) => this.standardGrav = standardGrav;

    public static OrbitPhysics StandardGravitationalParameter(OrbitGeometry geo, float standardGrav) => new OrbitPhysics(geo, standardGrav);

    public static OrbitPhysics CelestialMass(OrbitGeometry geo, float mass) => new OrbitPhysics(geo, mass * GRAV_CONST);

    public static OrbitPhysics OrbitalPeriod(OrbitGeometry geo, float orbitalPeriod)
    {
        //T=2pi*sqrt(a^3/u)
        //T^2 = 4*pi^2*a^3*u^-1
        //1 = u^-1 * 4*pi^2*a^3*T^-2
        //u = 4*pi^2*a^3*T^-2
        return new OrbitPhysics(geo, 4 * Square(Mathf.PI) * Cube(geo.SemiMajorAxis()) / Square(orbitalPeriod));
    }

    /// <summary>
    /// Standard Vis Viva Equation
    /// </summary>
    /// <param name="orbitalRadius"></param>
    /// <returns></returns>
    public float VelocityFromRadius(float orbitalRadius)
    {
        return Mathf.Sqrt(standardGrav * (2/orbitalRadius - 1/geo.SemiMajorAxis()));
    }


}
