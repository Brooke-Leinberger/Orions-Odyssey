using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System.Linq;

/// <summary>
/// Creates Geometry for orbiting bodies. Currently only works with closed orbits
/// TODO: add support for open orbits (issues with semi major axis being infinite or negative)
/// </summary>

public class OrbitGeometry
{

    #region Properties
    //-1 is used as a default value for most since almost every one should be non negative
    protected float

        //Identities (1)
        semiMajorAxis = -1,      //a
        eccentricity = -1,       //e

        //Principles (2)
        semiMinorAxis = -1,      //b
        linearEccentricity = -1, //c
        semiLatusRectum = -1,    //l
        latusScale = -1,         //s

        //Area as Time(5)
        fullArea = -1,
        scaleOfArea = -1,        //q
        dTheta = Mathf.Deg2Rad * 1f / 60f //resolution of lookup table (1 arcminute)
    ;

    enum Shape { CIRCLE, ELLIPSE, PARABOLA, HYPERBOLA }
    Shape shape;

    protected float[][] AnomalyAreaTable; //used to find the true anomaly, using a value equal to the area swept out by that true anomaly, from periapsis.
    private int polarity = -2;  //used when only the magnitude of the y coordinate certain, but polarity is important (may be deprecated) 
    #endregion

    #region Math Helpers
    private static float Square(float val) => val * val; //dont trust the pow function to be effecient enough
    private static float Cube(float val) => val * Square(val);
    private static float Mod(float val, float r) => (val % r + r) % r; //get true Modulo function
    #endregion

    public OrbitGeometry(float semiMajorAxis, float eccentricity)
    {
        Init(semiMajorAxis, eccentricity);
        if (eccentricity == 0) shape = Shape.CIRCLE;
        else if (0 < eccentricity && eccentricity < 1) shape = Shape.ELLIPSE;
        else if (eccentricity == 1) shape = Shape.PARABOLA;
        else if (eccentricity  > 1) shape = Shape.HYPERBOLA;
    }

    private void Init(float semiMajorAxis, float eccentricity)
    {
        //Indentites(1)
        this.semiMajorAxis = semiMajorAxis;
        this.eccentricity = eccentricity;
        PrincipleInit();
        AreaInit();
    }

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
        fullArea = semiMajorAxis * semiMinorAxis * Mathf.PI;
        scaleOfArea = semiMajorAxis * semiMinorAxis / 2; //The scale between a given area and its corresponding unit area
        dTheta = Mathf.Deg2Rad * (1f / 60f); //lookup table resolution is 1 arcsecond
        AnomalyAreaTableInit();
    }

    private void AnomalyAreaTableInit()
    {
        List<float[]> table = new List<float[]>();
        //iterate through the smallest resolution of true anomaly, then convert to eccentric to get unit area and store in lookup table.
        for (float i = 0; i <= 2 * Mathf.PI; i += dTheta) table.Add(new float[] { Area(i), i });
        AnomalyAreaTable = table.ToArray();
        //Can be taxing on machine in real time
        //Optimization concept: initialize only a small area ahead of the craft during high demand moments (e.g. changing orbit)
    }

    public float Area(float trueAnomaly) => 0.5f * semiMajorAxis * semiMinorAxis * UnitArea(trueAnomaly);

    public float UnitArea(float trueAnomaly)
    {
        float eccentricAnomaly = Mod(2 * Mathf.Atan(Mathf.Sqrt((1 - eccentricity) / (1 + eccentricity)) * Mathf.Tan(trueAnomaly / 2)), 2 * Mathf.PI);
        return eccentricAnomaly - eccentricity * Mathf.Sin(eccentricAnomaly);
    }

    public float LookUpTrueAnomalyTable(float area)
    { 
        area = Mod(area, fullArea);

        int indexUpperPot = AnomalyAreaTable.Length - 1, indexLowerPot = 0;
        int indexUpper, indexLower;
        bool inIndex;

        //binary search algorithm since it's already sorted (area function is continuously positive)
        do
        {
            indexUpper = indexUpperPot;
            indexLower = indexLowerPot;

            if (indexUpper - indexLower == 1) break;

            int middle = (indexLower + indexUpper) / 2;
            float median = AnomalyAreaTable[middle][0];
            if (area == median) return AnomalyAreaTable[middle][1];
            if (area < median) indexUpperPot = middle;
            else indexLowerPot = middle;

            //check if area will be out of bounds before commiting changes
            inIndex = AnomalyAreaTable[indexLowerPot][0] < area && area < AnomalyAreaTable[indexUpperPot][0];
        }
        while (inIndex);

        List<float[]> list = new List<float[]>();
        for (int i = indexLower; i < indexUpper + 1; i++) list.Add(AnomalyAreaTable[i]);
        return list.OrderBy(c => Mathf.Abs(area - c[0])).First()[1];
    }

    #region Accessors
    public float SemiMajorAxis() => semiMajorAxis;
    public float SemiMinorAxis() => semiMinorAxis;
    public float LinearEccentricity() => linearEccentricity;
    public float SemiLatusRectum() => semiLatusRectum;
    public float Eccentricity() => eccentricity;
    public float LatusScale() => latusScale;
    public float Area() => fullArea;
    public int Polarity() => polarity;

    #endregion

    #region Cartesian Conversions
    public float OrbitalRadius(float trueAnomaly)
    {
        return semiLatusRectum / (eccentricity * Mathf.Cos(trueAnomaly) + 1);
    }
    public float TrueAnomaly(float orbitalRadius)
    {
        return Mathf.Acos((semiLatusRectum - orbitalRadius)/(eccentricity * orbitalRadius));
    }
    public Vector2 CoordinatesAnomaly(float trueAnomaly)
    {
        return OrbitalRadius(trueAnomaly) * new Vector2(Mathf.Cos(trueAnomaly), Mathf.Sin(trueAnomaly));
    }
    public Vector2 CoordinatesRadius(float orbitalRadius)
    {
        return new Vector2(semiLatusRectum - orbitalRadius, Mathf.Sqrt(latusScale * (Square(linearEccentricity) - Square(semiMajorAxis - orbitalRadius)))) / eccentricity;
    }
    #endregion
}
