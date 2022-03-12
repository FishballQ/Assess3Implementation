using UnityEngine;

public class InWater : MonoBehaviour
{
    public bool customSurfaceLevel = true;
    public float surfaceLevel = 10f;

    private Collider coll;


    private void Start()
    {
        coll = GetComponent<Collider>();
    }

    //Get Surface level of the water
    public float GetYBound()
    {
        if (!customSurfaceLevel) 
            surfaceLevel = coll.bounds.max.y;
        return surfaceLevel;
    }
}