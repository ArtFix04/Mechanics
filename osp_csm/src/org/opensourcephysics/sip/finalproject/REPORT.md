# Project Summary: Packing and Floor-planning Simulation

**Name:** Artur Khachatryan
**Section:** B
**Student ID:** A09220183

This project implements and analyzes the dynamics of interacting particles using gravitational forces in a 2D environment. The simulation is developed in the **Open Source Physics (OSP)** Java-based environment and follows a series of tasks progressively improving performance and physical realism.

---

## Core Functionality

- Simulates both the **Packing Problem** (minimizing bounding box) and **Floor-planning Problem** (fixed rectangular area).
- Particles interact via **gravitational forces** using the **Verlet algorithm**, with **air resistance** added for numerical stability.
- Supports particle mass based on both **unit mass** and **unit density**.
- Includes **collision detection and resolution** to avoid overlap and simulate sticky behavior.

---

## Logging and Analysis

- Logs all simulation parameters, particle sizes, and both initial and final coordinates.
- Computes and stores total **potential energy**, **bounding box width/height**, and **bounding area** as a function of time.
- Saves results in structured text files (`log_initial.txt`, `log_final.txt`, `energy_area.txt`) for visualization and analysis.

---

## Optimization

- Improved computational efficiency using a **cell-based approximation** strategy.
- Interactions within the same cell are exact; distant cells are approximated by **super-particles** using their **center of mass** or **geometric center**.

---

## Outcome

The project provides a modular simulation framework for exploring spatial layout problems. It supports parameter studies, logs detailed outputs, and offers performance improvements.
