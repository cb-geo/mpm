using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System.IO;
using System.Linq;
using UnityEditor;
using System;
using Unity.VisualScripting;
using System.Diagnostics.Contracts;

[RequireComponent(typeof(MeshFilter))]
[RequireComponent(typeof(MeshRenderer))]
public class ParticleReader : MonoBehaviour
{
    // variable exposed to the editor to see the value during runtime
    public Vector3[] particles;
    public int frame = 0;
    public int frameCount = 1000;
    void Start()
    {
        TextAsset particlesInput = Resources.Load<TextAsset>("column3d/particles2");

        // creating unity mesh asset
        Mesh sceneMesh = new Mesh() { name = "Particles" };
        gameObject.GetComponent<MeshFilter>().mesh = sceneMesh;
        GetComponent<Renderer>().material.color = Color.red;

        string[] particleFileLines = particlesInput.text.Split("\n");

        // taking input from the first line of txt file
        int numParticles = int.Parse(particleFileLines[0].Split("\t")[0]);

        // putting the particle vertices from file to a Vec3 array
        particles = new Vector3[numParticles];
        int[] indices = new int[numParticles];

        for (int i = 1; i <= numParticles; i++)
        {
            particles[i - 1].Set(float.Parse(particleFileLines[i].Split("\t")[0]),
                              float.Parse(particleFileLines[i].Split("\t")[1]),
                              float.Parse(particleFileLines[i].Split("\t")[2]));
            indices[i - 1] = i - 1;
        }

        sceneMesh.vertices = particles;
        sceneMesh.SetIndices(indices, MeshTopology.Points, 0);
        sceneMesh.RecalculateBounds();

    }

    // Only works in scene view not in game view - Uncomment for debugging purpose
    private void OnDrawGizmos()
    {
        //if (particles == null)
        //    return;
        //for (int i = 0; i < particles.Length; i++)
        //{
        //    Gizmos.color = Color.blue;
        //    Gizmos.DrawSphere(gameObject.transform.position + particles[i], .01f);
        //}
    }

    private void FixedUpdate()
    {
        // load in the mpm generated geometry for the time steps of the simulation
        if (frame < frameCount)
        {
            int newNumParticles = System.IO.File.ReadLines("Assets/Resources/column3d/geometry_0" + frame.ToString("D3") + ".txt").Count() - 1;
            TextAsset updatedGeometry = Resources.Load<TextAsset>("column3d/geometry_0" + frame.ToString("D3"));
            string[] updatedGeometryFileLines = updatedGeometry.text.Split("\n");
            int[] indices = new int[newNumParticles];
            for (int i = 1; i <= newNumParticles; i++)
            {
                particles[i - 1].Set(float.Parse(updatedGeometryFileLines[i].Split(",")[0]),
                                  float.Parse(updatedGeometryFileLines[i].Split(",")[1]),
                                  float.Parse(updatedGeometryFileLines[i].Split(",")[2]));
                indices[i - 1] = i - 1;
            }
            gameObject.GetComponent<MeshFilter>().mesh.vertices = particles;
            gameObject.GetComponent<MeshFilter>().mesh.SetIndices(indices, MeshTopology.Points, 0);
            gameObject.GetComponent<MeshFilter>().mesh.RecalculateBounds();
            frame++;
        }
        else 
        { 

        }

    }

}
