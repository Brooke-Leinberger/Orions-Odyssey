using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class SatelliteMotion : MonoBehaviour
{
    Orbit satellite;
    public float SemiMajorAxis = 149.60e9f, Eccentricity = 0.0167086f, Mass = 1.989e30f, OrbitalPeriod = 31558149.7635456f; //defaults set to earth's orbit
    public float trueAnomaly = 0;
    public Orbit earth;
    public bool useOrbitalPeriod = false;

    // Start is called before the first frame update
    void Start()
    {
        if (!useOrbitalPeriod) earth = new Orbit(SemiMajorAxis, Eccentricity, Mass);
        else earth = Orbit.Period(SemiMajorAxis, Eccentricity, OrbitalPeriod);

        trueAnomaly = Mathf.Atan2(transform.localPosition.y, transform.localPosition.x);

    }

    // Update is called once per frame
    void Update()
    {
        float dTime = Time.deltaTime;
        trueAnomaly = earth.LookUpTrueAnomalyTable(trueAnomaly, dTime);
        Vector2 coor = earth.CoordinatesAnomaly(trueAnomaly);
        transform.position = new Vector3(coor.x, 0, coor.y);
        Debug.Log(1 / dTime);
    }

    
}
