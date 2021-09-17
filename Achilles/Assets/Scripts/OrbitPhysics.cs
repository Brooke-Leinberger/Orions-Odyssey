using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System.Linq;

public class OrbitPhysics : OrbitGeometry
{
    public const float GRAV_CONST = 6.6743e-11f; //physics
    private float
        standardGrav = -1,
        orbitalPeriod = -1,
        constantOfArea = -1,
        segmentedDTime = 1,
        segmentedDTheta = -1
        ;
 
    #region Math Helpers
    private static float Square(float val) => val * val; //dont trust the pow function to be effecient enough
    private static float Cube(float val) => val * Square(val);
    private static float Mod(float val, float r) => (val % r + r) % r; //get true Modulo function
    #endregion


    public OrbitPhysics(float semiMajorAxis, float eccentricity, float standardGrav, bool generateLookup = true) : base(semiMajorAxis, eccentricity, generateLookup)
    {
        this.standardGrav = standardGrav;
    }

    /// <summary>
    /// Standard Vis Viva Equation
    /// </summary>
    /// <param name="orbitalRadius"></param>
    /// <returns></returns>
    public float VelocityFromRadius(float orbitalRadius)
    {
        return Mathf.Sqrt(standardGrav * ( (2 / orbitalRadius) - (1 / semiMajorAxis) ) );
    }

    private static OrbitPhysics EmptyInit(OrbitGeometry geo)
    {
        return (OrbitPhysics)geo;
    }

    public OrbitPhysics StandardInit(float standardGrav)
    {
        //Area as Time (5)
        this.standardGrav = standardGrav;
        if(orbitalPeriod == -1) orbitalPeriod = 2 * Mathf.PI * Mathf.Sqrt(Cube(semiMajorAxis) / standardGrav); //time it takes to complete 1 revolution
        constantOfArea = Mathf.Sqrt(semiLatusRectum * standardGrav) / 2; //area swept out by satellite from orbital focus per second
        return this;
    }

    public OrbitPhysics MassInit(float celestialMass)
    {
        return StandardInit(celestialMass * GRAV_CONST);
    }

    public OrbitPhysics PeriodInit(float orbitalPeriod)
    {
        //T=2pi*sqrt(a^3/u)
        //T^2 = 4*pi^2*a^3*u^-1
        //1 = u^-1 * 4*pi^2*a^3*T^-2
        //u = 4*pi^2*a^3*T^-2
        this.orbitalPeriod = orbitalPeriod;
        StandardInit(4 * Square(Mathf.PI) * Cube(semiMajorAxis) / Square(orbitalPeriod));
        return this;
    }

    public OrbitPhysics SegmentedInit(float seconds)
    {
        if (!initComplete) AnomalyAreaTableInit(seconds * constantOfArea);
        return this;
    }

    public OrbitPhysics VelocityRadiusInit(float velocity, float orbitalRadius)
    {
        standardGrav = Square(velocity) / ((2 / orbitalRadius - 1) / semiMajorAxis);
        return this;
    }

    public float StandardGrav() => standardGrav;
    public float OrbitalPeriod() => orbitalPeriod;
    public float ConstantOfArea() => constantOfArea;

    public new OrbitPhysics SetStartAnomaly(float trueAnomaly)
    {
        startAnomaly = trueAnomaly;
        return this;
    }
}
