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

The `mesh` object defines the mesh (node and cells) and boundary conditions. The following inputs are considered. 

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
                    "dir": 1,
                    "sign_n": -1,
                    "friction": 0.5
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

The `particles` object defines material points from either a ASCII file or the built-in Gauss point generator. 

```JSON
{

    "particles": [
        {
            "generator": {
                "check_duplicates": true,
                "pset_id": 99,
                "material_id": 0,
                "particles_type": "P2D",
                "type": "file",
                "location": "particles.txt",
                "io_type": "Ascii2D",
            }
        },
        {
            "generator": {
                "check_duplicates": true,
                "pset_id": 98,
                "material_id": 1,
                "particle_type": "P2D",
                "type": "gauss",
                "nparticles_per_dir": 2, 
                "cset_id": -1
            }
        }
    ]

}
```

For both generators:
* `"check_duplicates"` flag is set `true` to check for repeating material points within the model.
* `"pset_id"` is the ID of the generated particle set.
* `"material_id"` corresponds to a specific material model defined below.
* `"particles_type"` is either `"P2D"` or `"P3D"`.

For particles from file:
* `"type"` is set to `"file"`.
* `"location"` defines the particle file name in the working directory. 
* `"io_type"` is either `"Ascii2D"` or `"Ascii3D"`.

For particles from Gauss generator:
* `"type"` is set to `"gauss"`.
* `"nparticles_per_dir"` is the number of material points per direction in either  
* `"cset_id"` defines which cells have material points generated; using `-1` generates material points in all cells. 

## Materials

The `materials` object defines that material types are used in analysis. Available materials include:
* Linear Elastic
* Mohr-Coulomb
* Modified Cam Clay
* Bonded NorSand

Check out material `.tcc` files for further details. 

## Material Sets

The `material_sets` object defines the relationship between `material_id` and `pset_id`. **This object is optional.** If provided, this relationship is used to redefine the `material_id` for each particle associated with the provided `pset_id`. This optional object allows for the particles from a single input file to be assigned multiple materials.

```JSON
{

    "materials_sets": [
        {
            "material_id": 0,
            "pset_id": 2
        }
    ]

}
```

## External Loading Conditions

The `external_loading_conditions` object specifies gravity, nodal forces, and particle forces. 

```JSON 
{

    "external_loading_conditions": {
        "gravity": [0.0, -9.81],
        "concentrated_nodal_forces": [
            {
                "nset_id": 0,
                "math_function_id": 0,
                "dir": 1,
                "force": 100
            }
        ],
        "particle_surface_traction": [
            {
                "pset_id": 1,
                "math_function_id": 1,
                "dir": 0,
                "traction": 10
            }
        ]
    }

}
```

For concentrated nodal forces:
* `"nset_id"` correspond to defined sets within the entity sets JSON file.
* `"math_function_id"` corresponds to a specific math function defined below.
* `"dir"` is the direction of the constraint `(0|1|2)`.
* `"force"` is the concentrated nodal force applied (consistent units).

For particle surface tractions:
* `"pset_id"` correspond to defined sets within the entity sets JSON file.
* `"math_function_id"` corresponds to a specific math function defined below.
* `"dir"` is the direction of the constraint `(0|1|2)`.
* `"traction"` is the surface traction applied (consistent units).

## Math Functions

The `math_function` object defines how loads change over time. 

```JSON
{

    "math_function": [
        {
            "id": 0,
            "type": "Linear",
            "xvalues": [0.0, 1.0, 10.0],
            "fxvalues": [0.0, 1.0, 1.0]
        }
    ]

}
```

* `"id"` corresponds the applied boundary loading above. 
* `"type"` is currently limited to linear functions. 
* `"xvalues"` is the model time. 
* `"fxvalues"` is the relative weight of the force or traction magnitude defined in the loading conditions. 

## Analysis

The `analysis` object defines the type of analysis, time step, and total number of steps.

```JSON
{

    "analysis": {
        "type": "MPMExplicit2D",
        "mpm_scheme": "usf",
        "velocity_update": true,
        "dt": 1.0e-5,
        "nsteps": 1001,
        "locate_particles": true,
        "nload_balance_steps": 100,
        "uuid": "output-name-000",
        "resume": {
            "resume": false,
            "uuid": "output-name-000",
            "step": 100
        },
        "damping": {
            "type": "Cundall",
            "damping_factor": 0.01
        }
    }

}
```

* `"velocity_update"` flag is set `true` to use nodal velocity to update particle velocity. 
* `"dt"` is the time step. 
* `"nsteps"` is the total number of steps. 
* `"locate_particles"` flag is set `false` to allow simulation to continue if material points leave the mesh domain. 
* `"nload_balance_steps"` is **optional** to redefine the dynamic load balancing; default is every 1000 steps. 
* `"uuid"` is the analysis name. 

### Type

Supported types:

| Analysis      | Description     |
|---------------|-----------------|
| MPMExplicit2D | explicit 2D MPM |
| MPMExplicit3D | explicit 3D MPM |

### Scheme

Support schemes:

| Schemes | Description                 |
|---------|-----------------------------|
| usf     | update stress first         |
| usl     | update stress last          |
| musl    | modified update stress last |

### Resume

Resume is optional.

### Damping

Damping is optional. Currently only Cundall damping can be specified. 

## Post Processing

The `post_processing` object defines what output files are generated. 

```JSON
{

    "post_processing":{
        "output_steps": 100,
        "path": "results/",
        "vtk": [ ],
        "vtk_statevars": [
            {
                "phase_id": 0,
                "statevars": [ ]
            }
        ]
    }

}
```

* `"output_steps"` is the number of steps between generated files. 
* `"path"` is where output files will be saved; path will be generated if it does not already exist. 

### VTK

Available VTK files include `stresses`, `strains`, `displacements`, and `velocities`.

### VTK State Variables

State variables are generated based on their phase ID (e.g., solid is `0`, liquid is `1`). Available state variables are material dependent. 


