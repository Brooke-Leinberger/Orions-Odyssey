using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class SatelliteMotion : MonoBehaviour
{
    Orbit satellite;
    public float SemiMajorAxis = 149.60e9f, Eccentricity = 0.0167086f, Mass = 1.989e30f, OrbitalPeriod = 31558149.7635456f; //defaults set to earth's orbit
    public float trueAnomaly = 0;
    public Orbit earth;
    public GameObject celestialBody = null;
    public bool useOrbitalPeriod = false;

    // Start is called before the first frame update
    void Start()
    {
        if (!useOrbitalPeriod) earth = new Orbit(SemiMajorAxis, Eccentricity, Mass);
        else earth = Orbit.Period(SemiMajorAxis, Eccentricity, OrbitalPeriod);

        trueAnomaly = Mathf.Atan2(transform.localPosition.y, transform.localPosition.x);
        Vector3 parent = celestialBody.transform.position;
        Vector3 center = new Vector3(parent.x - earth.LinearEccentricity(), 0, 0);
        Vector3 vertex = new Vector3(center.x - earth.SemiMajorAxis(), 0, 0);
        Vector3 covertex = new Vector3(center.x, 0, earth.SemiMinorAxis());
        Vector3 SemiLatus = new Vector3(parent.x, 0, earth.SemiLatusRectum());
        
        //draw debug lines for 1 hour
        Debug.DrawLine(parent, center, Color.magenta, 3600);
        Debug.DrawLine(center, vertex, Color.blue, 3600);
        Debug.DrawLine(center, covertex, Color.red, 3600);
        Debug.DrawLine(parent, SemiLatus, Color.green, 3600);

        Vector3 camera = center;
        camera.y = Camera.main.transform.position.y;
        Camera.main.transform.position = camera;

        for (int i = 0; i < 361; i++) Debug.Log(i +" : "+ earth.TrueToEccentric(i * Mathf.Deg2Rad) * Mathf.Rad2Deg);

    }

    // Update is called once per frame
    void Update()
    {
        float dTime = Time.deltaTime;
        trueAnomaly = earth.LookUpTrueAnomalyTable(trueAnomaly, dTime);
        Vector2 coor = earth.CoordinatesAnomaly(trueAnomaly);
        transform.position = new Vector3(coor.x, 0, coor.y);
    }

    
}
