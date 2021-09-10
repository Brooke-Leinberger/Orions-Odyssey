using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System.Linq;

public class OrbitGeometry : MonoBehaviour
{
    #region Properties
    //-1 is used as a default value for most since almost every one should be non negative
    private float

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
        orbitalPeriod = -1,      //T
        constantOfArea = -1,     //dA/dt
        scaleOfArea = -1,        //q
        dTheta = Mathf.Deg2Rad * 1f / 60f //resolution of lookup table
    ;

    enum Shape { CIRCLE, ELLIPSE, PARABOLA, HYPERBOLA }
    Shape shape;

    private float[][] AnomalyAreaTable; //used to find the true anomaly, using a value equal to the area swept out by that true anomaly, from periapsis.
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
        float eccentricAnomaly = Mod(2 * Mathf.Atan(Mathf.Sqrt((1 + eccentricity) / (1 - eccentricity)) * Mathf.Tan(trueAnomaly / 2)), 2 * Mathf.PI);
        return eccentricAnomaly - eccentricity * Mathf.Sin(eccentricAnomaly);
    }

    public float LookUpTrueAnomalyTable(float area)
    { 
        area = Mod(area, fullArea);
        float dAreaPlot = 2 * (fullArea / 2 - Area(Mathf.PI - dTheta / 2));
        //any two areas in the lookup table can't be farther away than the area change at apoapsis after (maxFrameRate ^ -1) seconds

        int indexUpperPot = AnomalyAreaTable.Length - 1, indexLowerPot = 0;
        int indexUpper, indexLower;
        bool inIndex;

        //binary search algorithm since it's already sorted (area function is continuously positive)
        do
        {
            indexUpper = indexUpperPot;
            indexLower = indexLowerPot;

            int middle = (indexLower + indexUpper) / 2;
            float median = AnomalyAreaTable[middle][0];
            if (area == median) return AnomalyAreaTable[middle][1];
            if (area < median) indexUpperPot = middle;
            else indexLowerPot = middle;

            //check if area will be out of bounds before commiting changes
            inIndex = AnomalyAreaTable[indexLowerPot][0] < area && area < AnomalyAreaTable[indexUpperPot][0];
        }
        while (inIndex);



        float smallestDifference = 2 * dAreaPlot;
        float[][] points = AnomalyAreaTable.Skip(indexLower).Take(indexUpper).ToArray();
        int index = -1;

        for (int i = 0; i < points.Length; i++)
        {
            float difference = Mathf.Abs(points[i][1] - area);
            if (Mathf.Abs(points[i][1] - area) < smallestDifference)
            {
                smallestDifference = difference;
                index = i;
            }

        }
        if (index == -1)
        {
            Debug.LogError("Could not find true anomaly corresponding to area: " + area);
            return -1;
        }
        else return AnomalyAreaTable[index][1];
    }

    #region Accessors
    public float SemiMajorAxis() => semiMajorAxis;
    public float SemiMinorAxis() => semiMinorAxis;
    public float LinearEccentricity() => linearEccentricity;
    public float SemiLatusRectum() => semiLatusRectum;
    public float Eccentricity() => eccentricity;
    public float LatusScale() => latusScale;
    public int Polarity() => polarity;
    public float ConstantOfArea() => constantOfArea;
    public float OrbitalPeriod() => orbitalPeriod;
    #endregion
}
