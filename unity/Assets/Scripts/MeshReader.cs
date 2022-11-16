using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System.IO;
using System.Linq;
using UnityEditor;
using System;
using Unity.VisualScripting;

[RequireComponent(typeof(MeshFilter))]
[RequireComponent(typeof(MeshRenderer))]
public class MeshReader : MonoBehaviour
{
    // variable exposed to the editor to see the value during runtime
    public Vector3[] mesh;

    void Start()
    {
        TextAsset meshInput = Resources.Load<TextAsset>("column3d/mesh");

        // creating unity mesh asset 
        Mesh sceneMesh = new Mesh() { name = "Mesh" };
        gameObject.GetComponent<MeshFilter>().mesh = sceneMesh;
        GetComponent<Renderer>().material.color = Color.green;

        string[] meshFileLines = meshInput.text.Split("\n");

        // taking input from the first line of txt file
        int numVertices = int.Parse(meshFileLines[0].Split("\t")[0]);

        // putting the mesh vertices from file to a Vec3 array
        mesh = new Vector3[numVertices];
        int[] indices = new int[numVertices];
        
        for (int i = 1; i <= numVertices; i++)
        {
            mesh[i - 1].Set(float.Parse(meshFileLines[i].Split("\t")[0]),
                              float.Parse(meshFileLines[i].Split("\t")[1]),
                              float.Parse(meshFileLines[i].Split("\t")[2]));
            indices[i - 1] = i - 1;
        }


        sceneMesh.vertices = mesh;
        sceneMesh.SetIndices(indices, MeshTopology.Points, 0);
        sceneMesh.RecalculateBounds();
    }


    // Only works in scene view not in game view - Uncomment for debugging purpose
    private void OnDrawGizmos()
    {
        //if (mesh == null)
        //    return;
        //for (int i = 0; i < mesh.Length; i++)
        //{
        //    Gizmos.color = Color.red;
        //    Gizmos.DrawSphere(gameObject.transform.position + mesh[i], .01f);
        //}
        //Gizmos.color = Color.yellow;
        //Gizmos.DrawWireCube(gameObject.transform.position + gameObject.GetComponent<MeshFilter>().sharedMesh.bounds.center, gameObject.GetComponent<MeshFilter>().sharedMesh.bounds.extents * 2);
    }

    private void Update()
    {

    }


}
