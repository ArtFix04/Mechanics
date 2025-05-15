package org.opensourcephysics.sip.finalproject;

import org.opensourcephysics.controls.*;
import org.opensourcephysics.frames.*;
import org.opensourcephysics.display.*;
import java.util.*;
import java.io.*;

public class PackingProblemOptimized extends AbstractSimulation {
    PlotFrame frame = new PlotFrame("x", "y", "Optimized Packing Problem");
    ArrayList<Particle> particles = new ArrayList<>();
    double dt;
    int N, T, step = 0;
    double drag = 0.01;
    double G = 1.0;
    double epsilonSq = 0.01;
    double cellSize = 5.0;
    Random rand = new Random();

    FileWriter initialLog, finalLog, energyLog;

    class Particle implements Drawable {
        double x, y, vx = 0, vy = 0, ax = 0, ay = 0, radius, mass;
        double ax_prev = 0, ay_prev = 0;
        Circle circle = new Circle();

        Particle(double x, double y, double radius) {
            this.x = x;
            this.y = y;
            this.radius = radius;
            this.mass = 1.0; // Use uniform mass or area * density
            circle.setXY(x, y);
            circle.pixRadius = (int) (radius * 10);
        }

        void saveAcceleration() {
            ax_prev = ax;
            ay_prev = ay;
        }

        void updatePosition() {
            x += vx * dt + 0.5 * ax_prev * dt * dt;
            y += vy * dt + 0.5 * ay_prev * dt * dt;
            circle.setXY(x, y);
        }

        void updateVelocity() {
            vx += 0.5 * (ax_prev + ax) * dt;
            vy += 0.5 * (ay_prev + ay) * dt;
            vx *= (1 - drag);
            vy *= (1 - drag);
        }

        public void draw(DrawingPanel panel, java.awt.Graphics g) {
            circle.draw(panel, g);
        }
    }

    class Cell {
        int row, col;
        Cell(int row, int col) {
            this.row = row;
            this.col = col;
        }
        public boolean equals(Object o) {
            if (!(o instanceof Cell)) return false;
            Cell c = (Cell) o;
            return row == c.row && col == c.col;
        }
        public int hashCode() {
            return Objects.hash(row, col);
        }
    }

    public void initialize() {
        frame.clearDrawables();
        particles.clear();
        step = 0;

        dt = control.getDouble("dt");
        N = control.getInt("N");
        T = control.getInt("T");
        drag = control.getDouble("drag");

        try {
            initialLog = new FileWriter("log_initial.txt");
            energyLog = new FileWriter("energy_area.txt");
            energyLog.write("Step\tPotentialEnergy\tWidth\tHeight\tArea\n");

            for (int i = 0; i < N; i++) {
                double radius = rand.nextDouble() * 3f;
                double x, y;
                boolean overlaps;
                do {
                    overlaps = false;
                    x = rand.nextDouble() * 20 + 2.5;
                    y = rand.nextDouble() * 20 + 2.5;
                    for (Particle p : particles) {
                        double dx = x - p.x;
                        double dy = y - p.y;
                        if (Math.sqrt(dx * dx + dy * dy) < radius + p.radius) {
                            overlaps = true;
                            break;
                        }
                    }
                } while (overlaps);

                Particle particle = new Particle(x, y, radius);
                particles.add(particle);
                frame.addDrawable(particle);
                initialLog.write(String.format("radius=%.2f x=%.4f y=%.4f\n", radius, x, y));
            }

            initialLog.write(String.format("Parameters: N=%d T=%d dt=%.4f drag=%.4f\n", N, T, dt, drag));
            initialLog.close();
        } catch (IOException e) {
            System.out.println("Initialization log failed: " + e.getMessage());
        }

        frame.setPreferredMinMax(0, 25, 0, 25);
    }

    void computeGravitationalAndStickyCollisions() {
        Map<Cell, List<Particle>> cellMap = new HashMap<>();

        for (Particle p : particles) {
            p.ax = 0;
            p.ay = 0;
            int row = (int) (p.y / cellSize);
            int col = (int) (p.x / cellSize);
            Cell c = new Cell(row, col);
            cellMap.computeIfAbsent(c, k -> new ArrayList<>()).add(p);
        }

        for (Map.Entry<Cell, List<Particle>> entryC : cellMap.entrySet()) {
            Cell cellC = entryC.getKey();
            List<Particle> particlesC = entryC.getValue();

            for (int i = 0; i < particlesC.size(); i++) {
                Particle pi = particlesC.get(i);
                for (int j = i + 1; j < particlesC.size(); j++) {
                    Particle pj = particlesC.get(j);
                    interactParticles(pi, pj);
                }
            }

            for (Map.Entry<Cell, List<Particle>> entryB : cellMap.entrySet()) {
                Cell cellB = entryB.getKey();
                if (cellB.equals(cellC)) continue;
                List<Particle> particlesB = entryB.getValue();
                double mx = 0, my = 0, totalMass = 0;
                for (Particle pb : particlesB) {
                    mx += pb.x * pb.mass;
                    my += pb.y * pb.mass;
                    totalMass += pb.mass;
                }

                double centerX, centerY;
                if (isAdjacent(cellC, cellB)) {
                    centerX = mx / totalMass;
                    centerY = my / totalMass;
                } else {
                    centerX = (cellB.col + 0.5) * cellSize;
                    centerY = (cellB.row + 0.5) * cellSize;
                }

                for (Particle pi : particlesC) {
                    applySuperParticleForce(pi, centerX, centerY, totalMass);
                }
            }
        }
    }

    boolean isAdjacent(Cell c1, Cell c2) {
        return Math.abs(c1.row - c2.row) <= 1 && Math.abs(c1.col - c2.col) <= 1;
    }

    void applySuperParticleForce(Particle p, double sx, double sy, double m) {
        double dx = sx - p.x;
        double dy = sy - p.y;
        double distSq = dx * dx + dy * dy + epsilonSq;
        double dist = Math.sqrt(distSq);
        double force = G * p.mass * m / distSq;
        double fx = force * dx / dist;
        double fy = force * dy / dist;
        p.ax += fx / p.mass;
        p.ay += fy / p.mass;
    }

    void interactParticles(Particle pi, Particle pj) {
        double dx = pj.x - pi.x;
        double dy = pj.y - pi.y;
        double distSq = dx * dx + dy * dy + epsilonSq;
        double dist = Math.sqrt(distSq);
        double minDist = pi.radius + pj.radius;
        double force = G * pi.mass * pj.mass / distSq;
        double fx = force * dx / dist;
        double fy = force * dy / dist;

        pi.ax += fx / pi.mass;
        pi.ay += fy / pi.mass;
        pj.ax -= fx / pj.mass;
        pj.ay -= fy / pj.mass;

        if (dist < minDist) {
            double overlap = minDist - dist;
            double correctionX = 0.5 * overlap * dx / dist;
            double correctionY = 0.5 * overlap * dy / dist;
            pi.x -= correctionX;
            pi.y -= correctionY;
            pj.x += correctionX;
            pj.y += correctionY;

            double totalVx = (pi.vx + pj.vx) / 2;
            double totalVy = (pi.vy + pj.vy) / 2;
            pi.vx = pj.vx = totalVx;
            pi.vy = pj.vy = totalVy;
        }
    }

    public void doStep() {
        if (step >= T) {
            stop();
            return;
        }

        for (Particle p : particles) {
            p.saveAcceleration();
        }

        for (Particle p : particles) {
            p.updatePosition();
        }

        computeGravitationalAndStickyCollisions();

        for (Particle p : particles) {
            p.updateVelocity();
        }

        logEnergyAndBounds();
        frame.repaint();
        step++;
    }

    void logEnergyAndBounds() {
        double energy = 0;
        double xmin = Double.POSITIVE_INFINITY, xmax = Double.NEGATIVE_INFINITY;
        double ymin = Double.POSITIVE_INFINITY, ymax = Double.NEGATIVE_INFINITY;

        for (int i = 0; i < particles.size(); i++) {
            Particle pi = particles.get(i);
            xmin = Math.min(xmin, pi.x - pi.radius);
            xmax = Math.max(xmax, pi.x + pi.radius);
            ymin = Math.min(ymin, pi.y - pi.radius);
            ymax = Math.max(ymax, pi.y + pi.radius);

            for (int j = i + 1; j < particles.size(); j++) {
                Particle pj = particles.get(j);
                double dx = pj.x - pi.x;
                double dy = pj.y - pi.y;
                double dist = Math.sqrt(dx * dx + dy * dy + epsilonSq);
                energy += -G * pi.mass * pj.mass / dist;
            }
        }

        double width = xmax - xmin;
        double height = ymax - ymin;
        double area = width * height;

        try {
            energyLog.write(String.format("%d\t%.6f\t%.4f\t%.4f\t%.4f\n", step, energy, width, height, area));
        } catch (IOException e) {
            System.out.println("Energy logging failed: " + e.getMessage());
        }
    }

    public void stop() {
        super.stop();
        try {
            finalLog = new FileWriter("log_final.txt");
            for (Particle p : particles) {
                finalLog.write(String.format("radius=%.2f x=%.4f y=%.4f\n", p.radius, p.x, p.y));
            }
            finalLog.close();
            energyLog.close();
            System.out.println("Simulation complete. Logs saved.");
        } catch (IOException e) {
            System.out.println("Final logging failed: " + e.getMessage());
        }
    }

    public void reset() {
        control.setValue("dt", 0.1);
        control.setValue("N", 10);
        control.setValue("T", 500);
        control.setValue("drag", 0.01);
        initialize();
    }

    public static void main(String[] args) {
        SimulationControl.createApp(new PackingProblemOptimized());
    }
}
