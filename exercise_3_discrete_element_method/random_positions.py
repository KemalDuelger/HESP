import random

LSYS = 3
radius = 0.05
N = 2000

positions = []
velocities = []

while len(positions) < N:
    x = random.uniform(radius, LSYS - radius)
    y = random.uniform(radius, LSYS - radius)
    z = random.uniform(radius, LSYS - radius)
    # Überlappung vermeiden:
    too_close = False
    for (px, py, pz) in positions:
        dist = ((x-px)**2 + (y-py)**2 + (z-pz)**2)**0.5
        if dist < 2*radius:
            too_close = True
            break
    if not too_close:
        positions.append((x, y, z))
        vx = random.uniform(-0.5, 0.5)
        vy = random.uniform(-0.5, 0.5)
        vz = random.uniform(-0.5, 0.5)
        velocities.append((vx, vy, vz))

with open("inputconfig_2000", "w") as f:
    f.write("num_time_step = 10000;\n")
    f.write("time_step_length = 0.0001;\n\n")
    f.write("gamma = 70;\n")
    f.write("K = 100000;\n")
    f.write(f"LSYS = {LSYS};\n")
    f.write(f"particle_num = {N};\n")
    f.write("cut_off_radius = 2.5;\n")
    f.write("gravity = -9.81;\n\n")

    f.write("particle_positions =  {\n")
    for p in positions:
        f.write(f"    ({p[0]:.4f},{p[1]:.4f},{p[2]:.4f}),\n")
    f.write("};\n\n")

    f.write("particle_velocities =  {\n")
    for v in velocities:
        f.write(f"    ({v[0]:.4f},{v[1]:.4f},{v[2]:.4f}),\n")
    f.write("};\n\n")

    f.write("particle_masses =  {\n")
    for _ in range(N):
        f.write("    1,\n")
    f.write("};\n\n")

    f.write("particle_radiuses =  {\n")
    for _ in range(N):
        f.write(f"    {radius},\n")
    f.write("};\n")