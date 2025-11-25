import os
import subprocess
import sys

import imageio.v3 as iio
import mujoco
import numpy as np

gammabot = """
<mujoco model="gammabot">
  <option integrator="RK4" cone="elliptic"/>
  <compiler inertiafromgeom="true" meshdir="assets"/>
  
  <default>
  <!-- friction = [sliding, torsional, rolling] -->
  <geom friction="0.005 0.0001 0.0001"
        condim="3"
        solref="0.01 1"
        solimp="0.99 0.999 0.001 0.5 2"/>
  </default>

  <asset>
    <texture name="grid" type="2d" builtin="checker" rgb1=".1 .2 .3"
     rgb2=".2 .3 .4" width="300" height="300"/>
    <material name="grid" texture="grid" texrepeat="8 8" reflectance=".2"/>
  </asset>
  
  <asset>
    <!-- gammabot.stl must be triangulated; MuJoCo assumes meters -->
    <mesh name="gammabot" file="gammabot.stl" scale="0.001 0.001 0.001"/>
  </asset>

  <worldbody>
    <geom size="1 1 .01" type="plane" material="grid" margin="0.001"/>
    <light pos="0 0 .6"/>
    <camera name="closeup" pos="0 -.3 .2" xyaxes="1 0 0 0 1 2"/>
    <body name="gammabot" pos="0 0 0">
      <!-- type=mesh uses the named mesh; add rgba to see it clearly -->
      <geom type="mesh" mesh="gammabot" density="1500" rgba="0.7 0.7 0.9 1" contype="1" 
      conaffinity="1"/>
      <freejoint/>
    </body>
  </worldbody>
  
  <keyframe>
    <key name="launch"
         qpos="-0.25 0 0  1 0 0 0"
         qvel="0 0 0   0 0 0"/>
  </keyframe>
  
  <visual>
    <global offwidth="1920" offheight="1088"/>
  </visual>
</mujoco>
"""

# noinspection PyArgumentList
model = mujoco.MjModel.from_xml_string(gammabot)
data = mujoco.MjData(model)

kid = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_KEY, "launch")
mujoco.mj_resetDataKeyframe(model, data, kid)

body_id = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_BODY, "gammabot")
F_body = np.array([0.005, 0.0, 0.0])  # force in body frame
T_body = np.array([0.0, 0.0, 0.0])  # torque in body frame
p_body = np.array([0.0, 0.0, 0.0])  # point of application wrt COM (body frame)
# # Pose
# data.qpos[0:3] = [0.0, 0.0, 10.70]  # x,y,z (m)
# data.qpos[3:7] = [1.0, 0.0, 0.0, 0.0]  # quaternion (w,x,y,z)
#
# # Velocity
# data.qvel[0:3] = [0.0, 0.4, 0.0]  # linear (m/s)
# data.qvel[3:6] = [0.0, 0.0, 2.0]  # angular (rad/s)

mujoco.mj_forward(model, data)
# with mujoco.Renderer(model) as renderer:
#     renderer.update_scene(data, camera="closeup")
#     img = renderer.render()
#
# iio.imwrite("gammabot.png", img)
#
# path = "gammabot.png"
# if sys.platform.startswith("win"):
#     os.startfile(path)
# elif sys.platform == "darwin":
#     subprocess.run(["open", path], check=False)
# else:
#     subprocess.run(["xdg-open", path], check=False)

duration = 3  # (seconds)
frame_rate = 60  # (Hz)

# Simulate and display video.
frames = []
mujoco.mj_resetDataKeyframe(model, data, 0)  # Reset the state to keyframe 0
with mujoco.Renderer(model, width=1920, height=1088) as renderer:
    while data.time < duration:
        if data.time >= 1.0:
            R_bw = data.xmat[body_id].reshape(3, 3)  # body->world rotation
            F_world = R_bw @ F_body
            # Equivalent torque about COM: rotate body torque + r × F (all in world)
            r_world = R_bw @ p_body
            T_world = (R_bw @ T_body) + np.cross(r_world, F_world)
            data.xfrc_applied[body_id][:] = np.hstack([F_world, T_world])
        else:
            data.xfrc_applied[body_id][:] = 0.0

        mujoco.mj_step(model, data)
        if len(frames) < data.time * frame_rate:
            renderer.update_scene(data, "closeup")
            pixels = renderer.render()
            frames.append(pixels)

iio.imwrite("gammabot.mp4", np.stack(frames, axis=0), fps=frame_rate)

if sys.platform.startswith("win"):
    os.startfile("gammabot.mp4")
elif sys.platform == "darwin":
    subprocess.run(["open", "gammabot.mp4"])
else:
    subprocess.run(["xdg-open", "gammabot.mp4"])
