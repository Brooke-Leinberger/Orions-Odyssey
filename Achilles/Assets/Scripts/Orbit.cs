using System.Collections;
using System.Collections.Generic;
using UnityEngine;

/// <summary>
/// 2D Orbit class for 2 body systems. Moves counter-clockwise.
/// Author: Brooke Leinberger
/// Date created: September 4th, 2021
/// </summary>
public class Orbit : MonoBehaviour
{
    public const float GRAV_CONST = 6.6743e-11f;
    enum orbitState {SUBORBITAL, ELLIPTICAL, PARABOLIC, HYPERBOLIC}
    private int polarity = 1; //used to force the position to either the top or bottom half of the orbit
    private float semiMajorAxis = -1, semiMinorAxis = -1, linearEccentricity = -1, semiLatusRectum = -1; //'a', 'b', 'c', and 'l' on documentation, respectively
    private float eccentricity = -1, latusScale = -1; //'S' on documentation
    private float orbitalPeriod = -1, constantOfArea = -1; //T, and dA/dt on documentation, respectively
    private float standardGrav = -1, planetaryRadius = -1;
    private float escapeVel = -1;

    #region Constructor Initializers
    private void Init(float semiMajorAxis, float eccentricity, float standardGrav)
    {
        this.semiMajorAxis = semiMajorAxis;
        this.eccentricity = eccentricity;
        this.standardGrav = standardGrav;
        latusScale = 1 - square(eccentricity);
        linearEccentricity = semiMajorAxis * eccentricity;
        orbitalPeriod = 2 * Mathf.PI * Mathf.Sqrt(semiMajorAxis);
        constantOfArea = Mathf.Sqrt(semiMajorAxis * semiLatusRectum * standardGrav) / 2;
    }

    private void InitMass(float semiMajorAxis, float eccentricity, float planetaryMass) => Init(semiMajorAxis, eccentricity, planetaryMass * GRAV_CONST);

    private void navInit(float orbitalRadius, float planetaryMass, float velocity, float velocityAnomaly)
    {
        semiMajorAxis = MajorVisViva(velocity, GRAV_CONST * planetaryMass, orbitalRadius);
        float eccen = (square(orbitalRadius) - (2 * orbitalRadius * semiMajorAxis)) * square(Mathf.Cos(velocityAnomaly - Mathf.PI / 2));
        eccen = (eccen / square(semiMajorAxis)) - 1;
        eccentricity = Mathf.Sqrt(eccen);
        //eccentricity = Mathf.Sqrt((square(orbitalRadius) - 2 * semiMajorAxis * orbitalRadius) * square(Mathf.Cos(velocityAnomaly - Mathf.PI/2))/square(semiMajorAxis) - 1);
        InitMass(semiMajorAxis, eccentricity, planetaryMass);
    }
    #endregion
    
    #region Constructors

    public Orbit(float semiMajorAxis, float eccentricity, float planetaryMass) => InitMass(semiMajorAxis, eccentricity, planetaryMass);

    public Orbit(float OrbitalRadius, float PlanetaryMass, float velocity, float velocityAnomaly) => navInit(OrbitalRadius, PlanetaryMass, velocity, velocityAnomaly);

    public Orbit(float orbitalRadius, float planetaryMass, Vector2 velocityVector)
    {
        float velocityAnomaly = Mathf.Atan2(velocityVector.y, velocityVector.x);
        navInit(orbitalRadius, planetaryMass, velocityVector.magnitude, velocityAnomaly);
    }
    #endregion


    #region Static Constructor functions
    public static float square(float val) => val * val; //dont trust the pow function to be effecient enough
    public static float cube(float val) => val * square(val);
    /// <summary>
    /// Uses a rearranged Vis-Viva Formula to get semi major axis from velocity, distance, and planetary mass
    /// </summary>
    /// <param name="velocity">The orbital speed of the satellite, in meters per second</param>
    /// <param name="mu">The standard gravitation parameter of the orbital focus</param>
    /// <param name="radius">The orbital radius</param>
    /// <returns></returns>
    public static float MajorVisViva(float velocity, float mu, float radius) => 1/(2 / radius - square(velocity) / mu);
    public static float VelocityVisViva(float semiMajorAxis, float mu, float radius) => Mathf.Sqrt(mu * (2 / radius - 1 / semiMajorAxis));
    public static float MuVisViva(float semiMajorAxis, float velocity, float radius) => square(velocity) / (2 / radius - 1 / semiMajorAxis);
    #endregion

    #region Navigational Formulae
    /// <summary>
    /// </summary>
    /// <param name="trueAnomaly">The angle the satellite has swept around the orbital focus from the periapsis, in degrees.</param>
    /// <returns>The distance from the orbital focus</returns>
    public float OrbitalRadius(float trueAnomaly) => semiLatusRectum / (eccentricity * Mathf.Cos(trueAnomaly) + 1);
    public float OrbitalRadiusDegrees(float trueAnomaly) => OrbitalRadius(Mathf.Deg2Rad * trueAnomaly);
    /// <summary>
    /// Returns true anomaly value between 0 and 2PI radians, velocity vector is needed for for more precision; corresponding possible negative true anomaly is either multiplied by -1, or subtracted from 360 degrees
    /// </summary>
    /// <param name="radius">The distance from the orbital focus; OrbitalRadius.</param>
    /// <returns>The angle the satellite has swept around the orbital focus from the periapsis, in radians.</returns>
    public float TrueAnomaly(float radius) => Mathf.Acos((semiLatusRectum - radius)/(radius * eccentricity));
    public float TrueAnomalyDegrees(float radius) => Mathf.Rad2Deg * TrueAnomaly(radius);

    /// <summary>
    /// </summary>
    /// <param name="radius">The distance from the orbital focus; OrbitalRadius.</param>
    /// <returns>The X Coordinate of the satellite in orbit, in the format shown in the radius derivation of the documentation</returns>
    public float XCoordinate(float radius) => (semiLatusRectum - radius) / eccentricity;
    /// <summary>
    /// 
    /// </summary>
    /// <param name="radius"></param>
    /// <returns></returns>
    public float AbsYCoordinate(float radius) => (latusScale * (square(linearEccentricity)-square(square(semiMajorAxis) - square(radius)))) / eccentricity;
    public float YCoordinate(float trueAnomaly) => AbsYCoordinate(OrbitalRadius(trueAnomaly)) * Mathf.Sign(Mathf.PI - trueAnomaly);
    public float YCoordinateDegrees(float trueAnomaly) => AbsYCoordinate(OrbitalRadiusDegrees(trueAnomaly)) * Mathf.Sign(Mathf.PI - trueAnomaly);
    public float Ycoordinate(float radius) => polarity * AbsYCoordinate(radius);
    //add velocity function as both magnitude-angle and vector2 data types.
    #endregion

    #region Accessors
    public float SemiMajorAxis() => semiMajorAxis;
    public float SemiMinorAxis() => semiMinorAxis;
    public float LinearEccentricity() => linearEccentricity;
    public float SemiLatusRectum() => semiLatusRectum;
    public float Eccentricity() => eccentricity;
    public float LatusScale() => latusScale;
    public int Polarity() => polarity;
    #endregion

    /// <summary>
    /// Sets polarity/sign of y coordinate in orbit, defaults to +1 if sign = zero.
    /// </summary>
    /// <param name="sign">determines if the y coordinate of satellite in the orbit class will be positive or negative. (Postive moving away, Negative moving back)</param>
    public void SetPolarity(int sign)
    {
        if (polarity != 0) polarity = System.Math.Sign(sign);
        else polarity = 1;
    }

}
