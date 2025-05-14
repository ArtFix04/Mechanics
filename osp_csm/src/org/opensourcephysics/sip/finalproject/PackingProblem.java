package org.opensourcephysics.sip.finalproject;

import org.opensourcephysics.controls.*;
import org.opensourcephysics.frames.*;
import org.opensourcephysics.display.*;
import java.util.*;

public class PackingProblem extends AbstractSimulation {
    PlotFrame frame = new PlotFrame("x", "y", "Gravitational Sticky Particles");
    ArrayList<Particle> particles = new ArrayList<>();
    double dt;
    int N;
    int T;
    int step = 0;
    double drag = 0.01;
    double G = 1.0;
    double epsilonSq = 0.01;
    Random rand = new Random();

    class Particle implements Drawable {
        double x, y, vx = 0, vy = 0, ax = 0, ay = 0, radius;
        double ax_prev = 0, ay_prev = 0;
        Circle circle = new Circle();

        Particle(double x, double y, double radius) {
            this.x = x;
            this.y = y;
            this.radius = radius;
            circle.setXY(x, y);
            circle.pixRadius = (int) (radius * 100);
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

    public void initialize() {
        frame.clearDrawables();
        particles.clear();
        step = 0;

        dt = control.getDouble("dt");
        N = control.getInt("N");
        T = control.getInt("T");
        drag = control.getDouble("drag");

        for (int i = 0; i < N; i++) {
            double radius = 0.1;
            double x, y;
            boolean overlaps;
            do {
                overlaps = false;
                x = rand.nextDouble() * 10;
                y = rand.nextDouble() * 10;
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
        }

        frame.setPreferredMinMax(0, 12, 0, 12);
    }

    void computeGravitationalAndStickyCollisions() {
        for (Particle p : particles) {
            p.ax = 0;
            p.ay = 0;
        }

        for (int i = 0; i < particles.size(); i++) {
            Particle pi = particles.get(i);
            for (int j = i + 1; j < particles.size(); j++) {
                Particle pj = particles.get(j);

                double dx = pj.x - pi.x;
                double dy = pj.y - pi.y;
                double distSq = dx * dx + dy * dy + epsilonSq;
                double dist = Math.sqrt(distSq);
                double minDist = pi.radius + pj.radius;

                // Gravitational force
                double force = G / distSq;
                double fx = force * dx / dist;
                double fy = force * dy / dist;

                pi.ax += fx;
                pi.ay += fy;
                pj.ax -= fx;
                pj.ay -= fy;

                // STICKY COLLISION HANDLING
                if (dist < minDist) {
                    // 1. Push particles apart to contact distance
                    double overlap = minDist - dist;
                    double correctionX = 0.5 * overlap * dx / dist;
                    double correctionY = 0.5 * overlap * dy / dist;
                    pi.x -= correctionX;
                    pi.y -= correctionY;
                    pj.x += correctionX;
                    pj.y += correctionY;

                    // 2. Set shared velocity (inelastic stick)
                    double totalVx = (pi.vx + pj.vx) / 2;
                    double totalVy = (pi.vy + pj.vy) / 2;
                    pi.vx = pj.vx = totalVx;
                    pi.vy = pj.vy = totalVy;
                }
            }
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

        frame.repaint();
        step++;
    }

    public void stop() {
        super.stop();
        System.out.println("Simulation complete.");
    }

    public void reset() {
        control.setValue("dt", 0.1);
        control.setValue("N", 10);
        control.setValue("T", 500);
        control.setValue("drag", 0.01);
        initialize();
    }

    public static void main(String[] args) {
        SimulationControl.createApp(new PackingProblem());
    }
}
