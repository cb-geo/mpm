// ========================================================================== //
// Copyright (c) 2014-2019 The University of Texas at Austin.                 //
// All rights reserved.                                                       //
//                                                                            //
// Licensed under the Apache License, Version 2.0 (the "License");            //
// you may not use this file except in compliance with the License.           //
// A copy of the License is included with this software in the file LICENSE.  //
// If your copy does not contain the License, you may obtain a copy of the    //
// License at:                                                                //
//                                                                            //
//     https://www.apache.org/licenses/LICENSE-2.0                            //
//                                                                            //
// Unless required by applicable law or agreed to in writing, software        //
// distributed under the License is distributed on an "AS IS" BASIS, WITHOUT  //
// WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.           //
// See the License for the specific language governing permissions and        //
// limitations under the License.                                             //
//                                                                            //
// ========================================================================== //

#include <iostream>

#include <cstdlib>
#include <math.h>
#include <string.h>
#include <sys/types.h>
#include <unistd.h>

#include "Renderer.h"
#include <Application.h>
#include <dtypes.h>

#include <ospray/ospray.h>

int mpiRank, mpiSize;

#include "Debug.h"

using namespace gxy;
using namespace std;

int samples_per_partition = 100;

#define WIDTH 1920
#define HEIGHT 1080

int width = WIDTH;
int height = HEIGHT;

float radius = 0.;

void syntax(char* a) {
  cerr << "syntax: " << a << " [options] data" << endl;
  cerr << "optons:" << endl;
  cerr << "  -D            run debugger" << endl;
  cerr << "  -n nsamples   number of samples in each partition (100)" << endl;
  cerr << "  -s x y        window size (" << WIDTH << "x" << HEIGHT << ")"
       << endl;
  cerr << "  -r radius     radius of samples (" << radius << ")" << endl;
  exit(1);
}

int run(int argc, char* argv[]) {
  string data = "";
  char* dbgarg;
  bool dbg = false;

  ospInit(&argc, (const char**)argv);

  Application theApplication(&argc, &argv);
  theApplication.Start();

  for (int i = 1; i < argc; i++) {
    if (argv[i][0] == '-') switch (argv[i][1]) {
        case 'D':
          dbg = true, dbgarg = argv[i] + 2;
          break;
        case 'n':
          samples_per_partition = atoi(argv[++i]);
          break;
        case 'r':
          radius = atof(argv[++i]);
          break;
        case 's':
          width = atoi(argv[++i]);
          height = atoi(argv[++i]);
          break;
        default:
          syntax(argv[0]);
      }
    else if (data == "")
      data = argv[i];
    else
      syntax(argv[0]);
  }

  Renderer::Initialize();
  theApplication.Run();

  RendererP theRenderer = Renderer::NewP();

  mpiRank = theApplication.GetRank();
  mpiSize = theApplication.GetSize();

  srand(mpiRank);

  Debug* d = dbg ? new Debug(argv[0], false, dbgarg) : NULL;

  // SampleMsg::Register();

  if (mpiRank == 0) {
    theRenderer->Commit();

    // create empty distributed container for volume data
    VolumeP volume = Volume::NewP();

    // import data to all processes, smartly distributes volume across
    // processses this import defines the partitioning of the data across
    // processses if subsequent Import commands have a different partition, an
    // error will be thrown

    volume->Import(data);
    volume->set_global_origin(-0.5, -0.5, -0.5);
    volume->set_local_offset(0.5, 0.5, 0.5);
    volume->set_deltas(0.5, 0.5, 0.5);
    volume->set_local_counts(6, 6, 6);

    int neighbors[6];
    neighbors[0] = -1;
    neighbors[1] = -1;
    neighbors[2] = -1;
    neighbors[3] = -1;
    neighbors[4] = -1;
    neighbors[5] = -1;

    float deltaX, deltaY, deltaZ;
    volume->get_deltas(deltaX, deltaY, deltaZ);
    std::cout << "Deltas: " << deltaX << "\t" << deltaY << "\t" << deltaZ
              << "\n";

    float ox, oy, oz;
    volume->get_local_origin(ox, oy, oz);
    std::cout << "Origin: " << ox << "\t" << oy << "\t" << oz << "\n";

    int nx, ny, nz;
    volume->get_local_counts(nx, ny, nz);
    std::cout << "Counts n: " << nx << "\t" << ny << "\t" << nz << "\n";

    // create empty distributed container for particles
    // particle partitioning will match volume partition
    ParticlesP samples = Particles::NewP();

    samples->CopyPartitioning(volume);
    // samples->allocate(samples_per_partition);

    Particle particle;
    for (int i = 0; i < 100; i++) {
      float x, y, z;

      x = ((float)rand() / RAND_MAX);
      y = ((float)rand() / RAND_MAX);
      z = ((float)rand() / RAND_MAX);
      particle.u.value = 1;
      particle.xyz.x = x;
      particle.xyz.y = y;
      particle.xyz.z = z;

      samples->push_back(particle);
    }

    /*
    // define action to perform on volume (see SampleMsg above)
    SampleMsg *smsg = new SampleMsg(volume, samples);
    smsg->Broadcast(true, true);
    */

    samples->Commit();

    DatasetsP theDatasets = Datasets::NewP();
    theDatasets->Insert("samples", samples);
    theDatasets->Commit();

    vector<CameraP> theCameras;

#if 0
    for (int i = 0; i < 20; i++)
    {
      CameraP cam = Camera::NewP();

      cam->set_viewup(0.0, 1.0, 0.0);
      cam->set_angle_of_view(45.0);

      float angle = 2*3.1415926*(i / 20.0);

      float vpx = 8.0 * cos(angle);
      float vpy = 8.0 * sin(angle);

      cam->set_viewpoint(vpx, vpy, 0.0);
      cam->set_viewdirection(-vpx, -vpy, 0.0);

      cam->Commit();
      theCameras.push_back(cam);
    }
#else
    CameraP cam = Camera::NewP();
    cam->set_viewup(0.0, 1.0, 0.0);
    cam->set_angle_of_view(45.0);
    cam->set_viewpoint(4.0, 0.0, 0.0);
    cam->set_viewdirection(-2.0, 0.0, 0.0);
    cam->Commit();
    theCameras.push_back(cam);
#endif

    ParticlesVisP pvis = ParticlesVis::NewP();
    pvis->SetName("samples");
    pvis->Commit(theDatasets);

    VisualizationP v = Visualization::NewP();
    v->AddVis(pvis);
    float light[] = {1.0, 2.0, 3.0};
    int t = 1;
    v->get_the_lights()->SetLights(1, light, &t);
    v->get_the_lights()->SetK(0.4, 0.6);
    v->get_the_lights()->SetShadowFlag(false);
    v->get_the_lights()->SetAO(0, 0.0);
    v->Commit(theDatasets);

    RenderingSetP theRenderingSet = RenderingSet::NewP();

    int indx = 0;
    for (auto c : theCameras) {
      RenderingP theRendering = Rendering::NewP();
      theRendering->SetTheOwner((indx++) % mpiSize);
      theRendering->SetTheSize(width, height);
      theRendering->SetTheCamera(c);
      theRendering->SetTheDatasets(theDatasets);
      theRendering->SetTheVisualization(v);
      theRendering->Commit();
      theRenderingSet->AddRendering(theRendering);
    }

    theRenderingSet->Commit();

    std::cerr << "RENDER\n";
    theRenderer->Render(theRenderingSet);
    // #ifdef GXY_WRITE_IMAGES
    std::cerr << "WAIT\n";
    theRenderingSet->WaitForDone();
    std::cerr << "WAIT DONE\n";

    // #endif

    theRenderingSet->SaveImages(string("samples"));

    theApplication.QuitApplication();
  }

  theApplication.Wait();
}
