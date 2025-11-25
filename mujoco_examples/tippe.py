import os
import subprocess
import sys

import imageio.v3 as iio
import mujoco
import numpy as np

tippe_top = """
<mujoco model="tippe top">
  <option integrator="RK4"/>

  <asset>
    <texture name="grid" type="2d" builtin="checker" rgb1=".1 .2 .3"
     rgb2=".2 .3 .4" width="300" height="300"/>
    <material name="grid" texture="grid" texrepeat="8 8" reflectance=".2"/>
  </asset>

  <worldbody>
    <geom size=".2 .2 .01" type="plane" material="grid"/>
    <light pos="0 0 .6"/>
    <camera name="closeup" pos="0 -.1 .07" xyaxes="1 0 0 0 1 2"/>
    <body name="top" pos="0 0 .02">
      <freejoint/>
      <geom name="ball" type="sphere" size=".02" />
      <geom name="stem" type="cylinder" pos="0 0 .02" size="0.004 .008"/>
      <geom name="ballast" type="box" size=".023 .023 0.005"  pos="0 0 -.015"
       contype="0" conaffinity="0" group="3"/>
    </body>
  </worldbody>

  <keyframe>
    <key name="spinning" qpos="0 0 0.02 1 0 0 0" qvel="0 0 0 0 1 200" />
  </keyframe>
</mujoco>
"""

# noinspection PyArgumentList
model = mujoco.MjModel.from_xml_string(tippe_top)
data = mujoco.MjData(model)

mujoco.mj_forward(model, data)
with mujoco.Renderer(model) as renderer:
    renderer.update_scene(data, camera="closeup")
    img = renderer.render()

iio.imwrite("tippe.png", img)

path = "tippe.png"  # or "tippe.jpg"
if sys.platform.startswith("win"):
    os.startfile(path)
elif sys.platform == "darwin":
    subprocess.run(["open", path], check=False)
else:
    subprocess.run(["xdg-open", path], check=False)

duration = 7  # (seconds)
frame_rate = 60  # (Hz)

# Simulate and display video.
frames = []
mujoco.mj_resetDataKeyframe(model, data, 0)  # Reset the state to keyframe 0
with mujoco.Renderer(model) as renderer:
    while data.time < duration:
        mujoco.mj_step(model, data)
        if len(frames) < data.time * frame_rate:
            renderer.update_scene(data, "closeup")
            pixels = renderer.render()
            frames.append(pixels)

iio.imwrite("tippe.mp4", np.stack(frames, axis=0), fps=frame_rate)

if sys.platform.startswith("win"):
    os.startfile("tippe.mp4")
elif sys.platform == "darwin":
    subprocess.run(["open", "tippe.mp4"])
else:
    subprocess.run(["xdg-open", "tippe.mp4"])
