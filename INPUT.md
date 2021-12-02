# Input JSON

## Overall Structure

Geomechanics MPM uses a `JSON` input file. The following file structure is used.

```JSON
{
    "title": "Exmaple Title",
    "mesh": {
        
    },
    "particles":[

    ],
    "materials":[

    ],
    "material_sets":[

    ],
    "external_loading_conditions":{

    },
    "math_functions": [

    ],
    "analysis": {

    },
    "post_processing": {

    }

}
```

## Mesh

The `mesh` object define the mesh (node and cells) and boundary conditions. The following inputs are considered. 

```JSON
{

    "mesh":{
        "mesh": "mesh.txt",
        "entity_sets": "entity_sets.json",
        "particles_volumes": "particles-volume.txt",
        "particles_stresses": "particles-stresses.txt",
        "particle_cells": "particles-cells.txt",
        "boundary_conditions":{
            "velocity_constraints": [
                {
                    "nset_id": 0,
                    "dir": 1,
                    "velocity": 0.0
                }
            ],
            "friction_constraints": [
                {
                    "nset_id": 1,
                    "dir" 1,
                    "sign_n": -1,
                    "friction" 0.5
                }
            ],
            "nodal_euler_angles": "nodal-euler-angles.txt",
            "particles_velocity_constraints": [
                {
                    "pset_id": 0,
                    "dir": 0,
                    "velocity": 0.0
                }
            ]
        },
        "cell_type": "ED2Q4",
        "isoparametric": false,
        "check_duplicates": false,
        "io_type": "Ascii3D",
        "node_type": "N2D"
    },

}
```

* `"isoparametric"` flag is set `true` for unstructured (non-prismatic) elementes. 
* `"check_diplicates"` flag is set `true` to check for repeating material points within the model. 
* `"io_type"` is either `"Ascii2D"` or `"Ascii3D"`.
* `"node_type"` is either `"N2D"` or `"N3D"`.

### Input Mesh Files

| File               | Description                                          | Optional? |
|--------------------|------------------------------------------------------|-----------|
| mesh               | nodal coordinates and cell connectivity              | no        |
| entity sets        | defines particle, node, and cell sets                | yes       |
| particles volumes  | defines particle volumes                             | yes       |
| particles stresses | defines initial particle stresses (Voigt convention) | yes       |
| particle cells     | defines initial guess of particle location           | yes       |


### Boundary Conditions

Optional boundary conditions include nodal friction constraints, velocity constraints, and Euler angles, as well as particle velocity constraints. 

For nodal and velocity constraints:
* `"nset_id"` and `"pset_id"` correspond to defined sets within the entity sets JSON file.
* `"dir"` is the direction of the constraint `(0|1|2)`.
* `"velocity"` is the specified velocity in declared direction.

For friction constraints:
* `"nest_id"` correspond to defined nodal sets within the entity sets JSON file.
* `"dir"` is the direction of resistance `(0|1|2)`.
* `"sign_n"` is the sign of the normal vector to the friction plane `(-1|+1)`.
* `"friction"` is the frictional coefficient. 

### Cell Type

The following cell types are currently supported. 

|Cell type  | Description                           |
|-----------|---------------------------------------|
|ED2T2      | 2D Triangle 3-noded element           |
|ED2T6      | 2D Triangle 3-noded element           |
|ED2Q4      | 2D Quadrilateral 4-noded element      |
|ED2Q8      | 2D Quadrilateral 8-noded element      |
|ED2Q9      | 2D Quadrilateral 9-noded element      |
|ED2Q16G    | 2D GIMP Quadrilateral 4-noded element |
|ED3H8      | 3D Hexahedron 8-noded element         |
|ED3H20     | 3D Hexahedron 20-noded element        |
|ED3H64     | 3D GIMP Hexahedron 64-noded element   |



## Particles

## Materials

## Material Sets

## External Loading Conditions

## Math Functions

## Analysis

## Post Processing

