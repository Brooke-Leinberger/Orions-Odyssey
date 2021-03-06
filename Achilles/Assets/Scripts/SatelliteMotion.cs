using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System.Linq;

public class SatelliteMotion : MonoBehaviour
{
    public float SemiMajorAxis = 149.60e9f, Eccentricity = 0.0167086f, Mass = 1.989e30f, OrbitalPeriod = 31558149.7635456f; //defaults set to earth's orbit
    public float trueAnomaly = 0;
    private OrbitPhysics earth;
    public GameObject celestialBody = null;
    public bool useOrbitalPeriod = false;

    public bool fullUpdate = false;
    public bool isCompleted;
    

    // Start is called before the first frame update
    void Start()
    {
        trueAnomaly = Mathf.Atan2(transform.localPosition.y, transform.localPosition.x);

        if (!useOrbitalPeriod) earth = new OrbitPhysics(SemiMajorAxis, Eccentricity, Mass);
        else earth = new OrbitPhysics(SemiMajorAxis, Eccentricity, -1, generateLookup: false)
                .SetStartAnomaly(trueAnomaly).PeriodInit(OrbitalPeriod);

        #region draw
        Vector3 parent = celestialBody.transform.position;
        Vector3 center = new Vector3(parent.x - earth.LinearEccentricity, 0, 0);
        Vector3 vertex = new Vector3(center.x - earth.SemiMajorAxis, 0, 0);
        Vector3 covertex = new Vector3(center.x, 0, earth.SemiMinorAxis);
        Vector3 SemiLatus = new Vector3(parent.x, 0, earth.SemiLatusRectum);

        //draw debug lines for 1 hour
        Debug.DrawLine(parent, center, Color.magenta, 3600);
        Debug.DrawLine(center, vertex, Color.blue, 3600);
        Debug.DrawLine(center, covertex, Color.red, 3600);
        Debug.DrawLine(parent, SemiLatus, Color.green, 3600);
        #endregion

        //center camera
        Vector3 camera = center;
        camera.y = Camera.main.transform.position.y;
        Camera.main.transform.position = camera;
    }

    // Update is called once per frame
    void Update()
    {

        if (fullUpdate) earth.AnomalyAreaTableInit(); //calculate whole trajectory
        else
        {
            earth.ClearLookupTable(); //clear work
            earth.SetStartAnomaly(trueAnomaly).SegmentedInit(1); //continually only calculate 1 second ahead of current position
        }

        isCompleted = earth.AreaTableComplete;
        Debug.Log(1 / Time.deltaTime);
        float area = earth.Area(trueAnomaly) + (earth.ConstantOfArea() * Time.deltaTime);
        trueAnomaly = earth.LookUpTrueAnomalyTable(area);
        Vector2 coor = earth.CoordinatesAnomaly(trueAnomaly);
        transform.position = new Vector3(coor.x, 0, coor.y);

    }


}
