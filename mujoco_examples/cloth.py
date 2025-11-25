import os
import subprocess
import sys

import imageio.v3 as iio
import mujoco
import numpy as np

cloth = """
<mujoco model="cloth_demo">
  <option timestep="0.002" gravity="0 0 -9.81" integrator="RK4"/>
  <visual>
    <quality shadowsize="4096"/>
  </visual>

  <default>
    <!-- Gentle damping helps cloth stability -->
    <joint damping="0.05"/>
    <geom friction="0.8 0.01 0.001" condim="6" solref="0.005 1" solimp="0.9 0.95 0.001"/>
  </default>

  <asset>
    <texture name="grid" type="2d" builtin="checker" rgb1=".15 .18 .22" rgb2=".25 .28 .32" width="512" height="512"/>
    <material name="matgrid" texture="grid" texrepeat="8 8" reflectance="0.1" rgba="1 1 1 1"/>
    <material name="matcarpet" rgba="0.85 0.85 0.9 1"/>
  </asset>

  <worldbody>
    <!-- Floor -->
    <geom name="floor" type="plane" size="5 5 0.05" material="matgrid" pos="0 0 0"/>

    <!-- Sphere to drape over -->
    <body name="ball" pos="0 0 0.35">
      <geom name="ball_geom" type="sphere" size="0.25" rgba="0.7 0.6 0.5 1"/>
    </body>

    <!-- Lights & camera -->
    <light pos="2 -2 3" dir="-1 1 -2"/>
    <camera name="iso" pos="0.0 -0.5 1.0" xyaxes="1 0 0  0 0.7 0.7"/>

    <!-- Cloth composite (your block) -->
    <flexcomp name="fabric" type="grid" dim="2" count="20 20 1" spacing="0.025 0.025 
    0.05" radius="0.01" pos="0 0 0.6">
      <!-- friction etc. if you want per-flex contact tuning -->
      <contact friction="8.0 0.1 0.01"/>
      <!-- pin two corners by grid coordinates (i j) -->
    </flexcomp>

  </worldbody>

  <!-- Let MuJoCo infer particle masses/inertias from geom size & density -->
  <compiler inertiafromgeom="true"/>
  
  <visual>
    <global offwidth="1920" offheight="1080"/>
  </visual>
</mujoco>
"""

# noinspection PyArgumentList
model = mujoco.MjModel.from_xml_string(cloth)
data = mujoco.MjData(model)

mujoco.mj_forward(model, data)

duration = 2  # (seconds)
frame_rate = 60  # (Hz)

# Simulate and display video.
frames = []
# mujoco.mj_resetDataKeyframe(model, data, 0)  # Reset the state to keyframe 0
with mujoco.Renderer(model, width=1920, height=1080) as renderer:
    next_t = 0.0
    dt_frame = 1.0 / frame_rate
    while data.time < duration:
        mujoco.mj_step(model, data)
        if data.time >= next_t:
            renderer.update_scene(data, "iso")
            frames.append(renderer.render())
            next_t += dt_frame

iio.imwrite("cloth.mp4", np.stack(frames, axis=0), fps=frame_rate)

if sys.platform.startswith("win"):
    os.startfile("cloth.mp4")
elif sys.platform == "darwin":
    subprocess.run(["open", "cloth.mp4"])
else:
    subprocess.run(["xdg-open", "cloth.mp4"])
