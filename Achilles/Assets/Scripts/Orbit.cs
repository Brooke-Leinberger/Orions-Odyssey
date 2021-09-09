using System.Collections;
using System.Collections.Generic;
using System.Linq;
using UnityEngine;

/// <summary>
/// 2D Orbit class for 2 body systems. Moves counter-clockwise.
/// </summary>
 
 /*
 Reference guide: Section numbers [§x.x.x(.x)] can used to see exact definitions, 
 and applications in the provided external documentation
 -first digit reference page, 
 -second digit references section, 
 -third digit refrerences either subsection or item number,
 -if there is a fourth digit, it refers to item number.
 
 (order of documentation subject to change; this is the order as of 9/9/2021)
 quick reference:
 
 Geometric Files:
 (1)  Preface and Axioms
 (2)  Principles Derivation
 (3)  Eccentricity Derivation //May change in placement
 (4)  Radius Derivation
 (5)  Area Derivation
 
 The Following Geometics Files are tentative/in progress:
 (6)  Parbolic and Hyperbolic Geometry
 (7)  True anomaly from Area Derivation
 (8)  Geometric Index
 
 The Following Physics Files are tentative/in progress
 Physics Files:
 (9)  Physics Derivation
 (10) Orbital functions as a function of time derivation
 (11) Astronaut Derivation of Properties
 (12) Changes in orbital characteristics
 (13) Patched Conics Derivation and Application
 (14) Physics Index
 (15) Programming Index //almost all functions in this fill will reference this
 (16) Complete Index
 
 Author: Brooke Leinberger
 Date created: September 4th, 2021
 */
public class Orbit : MonoBehaviour
{
    #region Properties
    public const float GRAV_CONST = 6.6743e-11f; //physics
    private int polarity = 1;
    /* 
     * There are many ways to construct an orbit, the documentation can show many many ways to convert values to
     * to useful parameters that can reconstruct the geometry and/or the physics of an orbit.
     * 
     * For this class, we will have our base constructor use :
     * - Semi Major Axis and Eccentricity for geometry (§1.2.4 & §1.2.7, respectively)
     * - Standard Gravitation Parameter for Physics (§[WIP])
     * 
     * Constructors accepting other parameters will first convert those parameters to the parameters above before 
     * reconstructing the orbit.
     */

    //Negative one is used as a default value for most properties since it should be impossible to naturally get.
    private float
        //Identities (1)
        semiMajorAxis = -1,      //a
        eccentricity = -1,       //e

        //Principles (2)
        semiMinorAxis = -1,      //b
        linearEccentricity = -1, //c
        semiLatusRectum = -1,    //l
        latusScale = -1,         //s

        //Eccentricity Derivation (3)

        //Area as Time(5)
        orbitalPeriod = -1,      //T
        constantOfArea = -1,     //dA/dt
        scaleOfArea = -1,       //q

        //Physics (TBD)
        standardGrav = -1; //physics

    private float maxFrameRate = 60; //given as frames per second, used to calculate resolution of lookup table.
    private float dTheta = -1; //measures the resolution of the Area_Anomaly lookup table
    private float[][] areaAnomalyPlot = null; //lookup table

    #endregion

    #region Constructor Initializers

    private void PrincipleInit()
    {
        //Principles(2)
        latusScale = 1 - Square(eccentricity);
        linearEccentricity = semiMajorAxis * eccentricity;
        semiMinorAxis = semiMajorAxis * Mathf.Sqrt(latusScale);
        semiLatusRectum = semiMajorAxis * latusScale;
    }
    private void AreaInit()
    {
        //Area as Time (5)
        orbitalPeriod = 2 * Mathf.PI * Mathf.Sqrt(Cube(semiMajorAxis) / standardGrav); //time it takes to complete 1 revolution
        constantOfArea = Mathf.Sqrt(semiLatusRectum * standardGrav) / 2; //area swept out by satellite from orbital focus per second
        scaleOfArea = semiMajorAxis * semiMinorAxis / 2; //The scale between a given area and its corresponding unit area
        dTheta = Mathf.Deg2Rad * (1f / 60f); //lookup table resolution is 1 arcsecond
        LookUpTableInit();
    }
    private void LookUpTableInit()
    {
        List<float[]> table = new List<float[]>();
        //iterate through the smallest resolution of true anomaly, then convert to eccentric to get unit area and store in lookup table.
        for (float i = 0; i <= 2 * Mathf.PI; i += dTheta) table.Add(new float[] { ScaledArea(TrueToEccentric(i)), i });
        areaAnomalyPlot = table.ToArray();
    }

    //Every constructor eventually references Init
    /// <summary>
    /// Initializes the standard orbit object using the most basic identites.
    /// All other initializing functions are based off this one.
    /// </summary>
    /// <param name="semiMajorAxis"></param>
    /// <param name="eccentricity"></param>
    /// <param name="standardGrav"></param>
    private void Init(float semiMajorAxis, float eccentricity, float standardGrav)
    {
        //Indentites(1)
        this.semiMajorAxis = semiMajorAxis;
        this.eccentricity = eccentricity;
        this.standardGrav = standardGrav;
        PrincipleInit();
        AreaInit();
    }
    /// <summary>
    /// Substitutes Standard Gravitation Parameter for Planetary Mass in the standard init function
    /// </summary>
    /// <param name="semiMajorAxis"></param>
    /// <param name="eccentricity"></param>
    /// <param name="planetaryMass"></param>
    private void InitMass(float semiMajorAxis, float eccentricity, float planetaryMass)
    {
        //Standard Gravitational parameter is the product of planetary mass and the gravitational constant.
        Init(semiMajorAxis, eccentricity, planetaryMass * GRAV_CONST);
    }
    /// <summary>
    /// Substitutes Standard Gravitational Period for Orbital Period and Semi Major Axis
    /// </summary>
    /// <param name="semiMajorAxis"></param>
    /// <param name="eccentricity"></param>
    /// <param name="orbitalPeriod"></param>
    private void InitPeriod(float semiMajorAxis, float eccentricity, float orbitalPeriod)
    {
        //rearrange Kepler's 3rd Law for orbital period for mu to get standard gravitational parameter
        Init(semiMajorAxis, eccentricity, MuOrbitalPeriod(semiMajorAxis, orbitalPeriod));
    }
    /// <summary>
    /// Derives Orbital Identites through navigational data
    /// </summary>
    /// <param name="orbitalRadius"></param>
    /// <param name="planetaryMass"></param>
    /// <param name="velocity"></param>
    /// <param name="velocityAnomaly"></param>
    private void NavInt(float orbitalRadius, float planetaryMass, float velocity, float velocityAnomaly)
    {
        semiMajorAxis = MajorVisViva(velocity, GRAV_CONST * planetaryMass, orbitalRadius);
        float eccen = (Square(orbitalRadius) - (2 * orbitalRadius * semiMajorAxis)) * Square(Mathf.Cos(velocityAnomaly - Mathf.PI / 2));
        eccen = (eccen / Square(semiMajorAxis)) - 1;
        eccentricity = Mathf.Sqrt(eccen);
        //eccentricity = Mathf.Sqrt((square(orbitalRadius) - 2 * semiMajorAxis * orbitalRadius) * square(Mathf.Cos(velocityAnomaly - Mathf.PI/2))/square(semiMajorAxis) - 1);
        InitMass(semiMajorAxis, eccentricity, planetaryMass);
    }
    #endregion

    #region Constructors

    private Orbit()
    {
    }

    public Orbit (float semiMajorAxis, float eccentricity, float planetaryMass) => InitMass(semiMajorAxis, eccentricity, planetaryMass);

    public Orbit(float OrbitalRadius, float PlanetaryMass, float velocity, float velocityAnomaly) => NavInt(OrbitalRadius, PlanetaryMass, velocity, velocityAnomaly);

    public Orbit (float orbitalRadius, float planetaryMass, Vector2 velocityVector)
    {
        float velocityAnomaly = Mathf.Atan2(velocityVector.y, velocityVector.x);
        NavInt(orbitalRadius, planetaryMass, velocityVector.magnitude, velocityAnomaly);
    }
    public static Orbit Period(float semiMajorAxis, float eccentricity, float orbitalPeriod)  
    {
        Orbit init = new Orbit();
        init.InitPeriod(semiMajorAxis, eccentricity, orbitalPeriod);
        return init;
    }
    #endregion

    #region Static Constructor functions
    private static float Square(float val) => val * val; //dont trust the pow function to be effecient enough
    private static float Cube(float val) => val * Square(val);
    private static float Mod(float val, float r) => (val % r + r) % r; //get true Modulo function


    //Phyisics Functions:
    #region Vis Viva Functions
    //The following functions with "VisViva" in the name are all rearanged formulations of the Vis-Viva formula
    /// <summary>
    /// Used to derive Speed (no direction information directly included) from various observational information
    /// </summary>
    /// <param name="semiMajorAxis">; §1.2.4</param>
    /// <param name="mu"></param>
    /// <param name="radius"></param>
    /// <returns></returns>
    public static float VelocityVisViva(float semiMajorAxis, float mu, float radius)
    {
        return Mathf.Sqrt(mu * ((2 / radius) - 1 / (semiMajorAxis)));
    }
    /// <summary>
    /// Uses a rearranged Vis-Viva Formula to get semi major axis from velocity, distance, and planetary mass
    /// </summary>
    /// <param name="velocity">The orbital speed of the satellite, in meters per second</param>
    /// <param name="mu">The standard gravitation parameter of the orbital focus</param>
    /// <param name="radius">The orbital radius</param>
    /// <returns></returns>
    public static float MajorVisViva(float velocity, float mu, float radius)
    {
        return 1 / ((2 / radius) - (Square(velocity) / mu));
    }
    public static float MuVisViva(float semiMajorAxis, float velocity, float radius) => Square(velocity) / ((2 / radius - 1) / semiMajorAxis);
    public static float MuOrbitalPeriod(float semiMajorAxis, float orbitalPeriod) => 4 * Square(Mathf.PI) * Cube(semiMajorAxis) / Square(orbitalPeriod);
    #endregion

    #endregion

    #region Navigational Formulae
    /// <summary>
    /// </summary>
    /// <param name="trueAnomaly">The angle the satellite has swept around the orbital focus from the periapsis, in degrees.</param>
    /// <returns>The distance from the orbital focus</returns>
    public float OrbitalRadius(float trueAnomaly) => semiLatusRectum / (eccentricity * Mathf.Cos(trueAnomaly) + 1);
    /// <summary>
    /// Returns true anomaly value between 0 and 2PI radians, velocity vector is needed for for more precision; corresponding possible negative true anomaly is either multiplied by -1, or subtracted from 360 degrees
    /// </summary>
    /// <param name="radius">The distance from the orbital focus; OrbitalRadius.</param>
    /// <returns>The angle the satellite has swept around the orbital focus from the periapsis, in radians.</returns>
    public float TrueAnomaly(float radius) => Mathf.Acos((semiLatusRectum - radius)/(radius * eccentricity));


    /// <summary>
    /// </summary>
    /// <param name="radius">The distance from the orbital focus; OrbitalRadius.</param>
    /// <returns>The X Coordinate of the satellite in orbit, in the format shown in the radius derivation of the documentation</returns>
    public float XCoordinate(float orbitalRadius) => (semiLatusRectum - orbitalRadius) / eccentricity;
    /// <summary>
    /// 
    /// </summary>
    /// <param name="radius"></param>
    /// <returns></returns>
    public float AbsYCoordinate(float orbitalRadius) => (latusScale * (Square(linearEccentricity)-Square(Square(semiMajorAxis) - Square(orbitalRadius)))) / eccentricity;
    /// <summary>
    /// Warning: Make sure polarity property is updated before calling
    /// </summary>
    /// <param name="orbitalRadius"></param>
    /// <returns></returns>
    public float YCoordinate(float orbitalRadius) => polarity * AbsYCoordinate(orbitalRadius);

    public float XCoordinateAnomaly(float trueAnomaly) => OrbitalRadius(trueAnomaly) * Mathf.Cos(trueAnomaly);
    public float YCoordinateAnomaly(float trueAnomaly) => OrbitalRadius(trueAnomaly) * Mathf.Sin(trueAnomaly);

    public Vector2 AbsCoordinates(float orbitalRadius) => new Vector2(XCoordinate(orbitalRadius), AbsYCoordinate(orbitalRadius));
    public Vector2 CoordinatesAnomaly(float trueAnomaly) => new Vector2(XCoordinateAnomaly(trueAnomaly), YCoordinateAnomaly(trueAnomaly));
    /// <summary>
    /// Warning: Make sure polarity property is updated before calling
    /// </summary>
    /// <param name="orbitalRadius"></param>
    /// <returns></returns>
    public Vector2 Coordinates(float orbitalRadius) => new Vector2(XCoordinate(orbitalRadius), YCoordinate(orbitalRadius));
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
    public float StandardGrav() => standardGrav;
    public float ConstantOfArea() => constantOfArea;
    public float OrbitalPeriod() => orbitalPeriod;
    #endregion

    #region Sets
    /// <summary>
    /// Sets polarity/sign of y coordinate in orbit, defaults to +1 if sign = zero.
    /// </summary>
    /// <param name="sign">determines if the y coordinate of satellite in the orbit class will be positive or negative. (Postive moving away, Negative moving back)</param>
    public void SetPolarity(int sign)
    {
        if (polarity != 0) polarity = System.Math.Sign(sign);
        else polarity = 1;
    }
    #endregion

    #region Areas
    //Area as Time (5):
    //Conversions between Eccentric Anomaly and True Anomaly
    public float TrueToEccentric(float trueAnomaly) => Mod(2 * Mathf.Atan(Mathf.Sqrt((1 - eccentricity) / (1 + eccentricity)) * Mathf.Tan(trueAnomaly / 2)), 2 * Mathf.PI);
    public float EccentricToTrue(float eccentricAnomaly) => Mod(2 * Mathf.Atan(Mathf.Sqrt((1 + eccentricity) / (1 - eccentricity)) * Mathf.Tan(eccentricAnomaly / 2)), 2 * Mathf.PI);
    //Determine Area from Eccentric Anomaly
    /// <summary>
    /// "f(theta)" in the Documentation, (5) Area as Time
    /// </summary>
    /// <param name="eccentricAnomaly"></param>
    /// <returns></returns>
    public float UnitArea(float eccentricAnomaly)
    {
        float area = eccentricAnomaly - eccentricity * Mathf.Sin(eccentricAnomaly);
        return area;
    }
    public float ScaledArea(float eccentricAnomaly) => 0.5f * semiMajorAxis * semiMinorAxis * UnitArea(eccentricAnomaly);
    public float DeltaTheta() =>  Mathf.Atan((VelocityVisViva(semiMajorAxis, standardGrav, semiMajorAxis + linearEccentricity) / maxFrameRate) / (semiMajorAxis + linearEccentricity));
    public void DeltaTheta(float refreshRate) =>  Mathf.Atan((VelocityVisViva(semiMajorAxis, standardGrav, semiMajorAxis + linearEccentricity) / refreshRate) / (semiMajorAxis + linearEccentricity));
    public float fullArea() => semiMajorAxis * semiMinorAxis * Mathf.PI;
    public float LookUpTrueAnomalyTable(float theta1, float deltaTime)
    {
        float dAreaPlot = 2 * (fullArea()/2 - ScaledArea(Mathf.PI - dTheta / 2));
        float deltaArea = deltaTime * constantOfArea;
        float area1 = ScaledArea(TrueToEccentric(theta1));
        float area2 = Mod(area1 + deltaArea, fullArea());
        //any two areas in the lookup table can't be farther away than the area change at apoapsis after (maxFrameRate ^ -1) seconds\
        float[][] points = areaAnomalyPlot.Where(c => Mathf.Abs(c[0] - area2) < dAreaPlot * 2).ToArray();

        float smallestDifference = 2 * dAreaPlot;
        int index = 0;


        //Debug.Log("Points.Count: " + points.Count());
        //Debug.Log(ScaledArea(181 * Mathf.Deg2Rad));

        for (int i = 0; i < points.Length; i++)
        {
            //Debug.Log(Mathf.Rad2Deg * points[i][1] + " : " + points[i][0] + " : " + Mathf.Rad2Deg * EccentricToTrue(points[i][1]));
            float difference = Mathf.Abs(points[i][1] - area2);
            if (Mathf.Abs(points[i][1] - area2) < smallestDifference)
            {
                smallestDifference = difference;
                index = i;
            }
            //Debug.Log(Mathf.Rad2Deg * points[i][1]);
            
        }
        return points[index][1];
    }
    #endregion

    public string Log()
    {
        string log = "";
        log += "Semi Major Axis: " + semiMajorAxis + "\n";
        log += "Eccentricity: " + eccentricity + "\n";
        log += "Standard Gravitational Parameter: " + standardGrav + "\n";
        log += "Orbital Period: " + orbitalPeriod + "\n";
        log += "Area covered / second: " + constantOfArea + "\n";
        log += "dTheta: " + dTheta + "\n";
        return log;
    }
}
