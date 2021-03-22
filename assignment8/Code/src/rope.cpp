#include <iostream>
#include <vector>

#include "CGL/vector2D.h"

#include "mass.h"
#include "rope.h"
#include "spring.h"

namespace CGL {

    Rope::Rope(Vector2D start, Vector2D end, int num_nodes, float node_mass, float k, vector<int> pinned_nodes)
    {
        // TODO (Part 1): Create a rope starting at `start`, ending at `end`, and containing `num_nodes` nodes.

        
//        Comment-in this part when you implement the constructor
        Vector2D step = (start - end) / (num_nodes - 1);
        for(int i = 0; i != num_nodes ; ++i)
        {
            Vector2D pos = start + step * i;
            masses.push_back(new Mass (pos, node_mass, false));
            if(i >= 1)
            springs.push_back(new Spring(masses[i-1], masses.back(), k));
        }
        for (auto &i : pinned_nodes) 
        {
            masses[i]->pinned = true;
        }
    }

    void Rope::simulateEuler(float delta_t, Vector2D gravity)
    {
        // std::vector<Vector2D> masses; 
        for (auto &s : springs)
        {
            // TODO (Part 2): Use Hooke's law to calculate the force on a node
            Vector2D a_b = s->m1->position - s->m2->position;
            float dis = a_b.norm();
            a_b /= dis;
            Vector2D fab =  - s->k * (dis -  s->rest_length) * a_b ;
            s->m1->forces += fab;
            s->m2->forces += - fab;
        }
        for (auto &m : masses)
        {
            if (!m->pinned)
            {
                // TODO (Part 2): Add the force due to gravity, then compute the new velocity and position
                // acceleration
                
                Vector2D ac = (m->forces + gravity - 0.01 * m->velocity) / m->mass;
                Vector2D v_next = m->velocity + ac * delta_t;
                
                // explicit
                // m->position += m->velocity * delta_t;
                // implicite
                m->position += v_next * delta_t;
                m->velocity= v_next;

                // TODO (Part 2): Add global damping

            }

            // Reset all forces on each mass
            m->forces = Vector2D(0, 0);
        }

    }
    void Rope::simulateEulerExplicit(float delta_t, Vector2D gravity)
    {
        // std::vector<Vector2D> masses;
        delta_t /= 5; 
        for (auto &s : springs)
        {
            // TODO (Part 2): Use Hooke's law to calculate the force on a node
            Vector2D a_b = s->m1->position - s->m2->position;
            float dis = a_b.norm();
            a_b /= dis;
            Vector2D fab =  - s->k * (dis -  s->rest_length) * a_b ;
            s->m1->forces += fab;
            s->m2->forces += - fab;
        }
        for (auto &m : masses)
        {
            if (!m->pinned)
            {
                // TODO (Part 2): Add the force due to gravity, then compute the new velocity and position
                // acceleration
                
                Vector2D ac = (m->forces + gravity - 0.0005 * m->velocity) / m->mass;
                Vector2D v_next = m->velocity + ac * delta_t;
                
                // explicit
                m->position += m->velocity * delta_t;
                // implicite
                // m->position += v_next * delta_t;
                m->velocity= v_next;

                // TODO (Part 2): Add global damping

            }

            // Reset all forces on each mass
            m->forces = Vector2D(0, 0);
        }

    }
    void Rope::simulateVerlet(float delta_t, Vector2D gravity)
    {
        for (auto &s : springs)
        {
            // TODO (Part 3): Simulate one timestep of the rope using explicit Verlet （solving constraints)
            Vector2D a_b = s->m1->position - s->m2->position;
            float dis = a_b.norm();
            a_b /= dis;
            Vector2D fab =  - s->k * (dis -  s->rest_length) * a_b ;
            s->m1->forces += fab;
            s->m2->forces += - fab;

        }

        for (auto &m : masses)
        {
            if (!m->pinned)
            {

                // m->position += v_next * delta_t;
                // m->velocity= v_next;
                Vector2D tmp_position = m->position;
                Vector2D ac = (m->forces + gravity) / m->mass;
                m->position += (1 - 0.00005) * (m->position - m->start_position) + ac * delta_t * delta_t;
                m->start_position = tmp_position;
                m->forces = Vector2D(0, 0);
                // TODO (Part 4): Add global Verlet damping
                
            }
        }
    }
    void Rope::simulateVerletConstraints(float delta_t, Vector2D gravity)
    {
        for (auto &s : springs)
        {
            // TODO (Part 3): Simulate one timestep of the rope using explicit Verlet （solving constraints)
            Vector2D a_b = s->m1->position - s->m2->position;
            float dis = a_b.norm();
            a_b /= dis;
            Vector2D constra = - 0.5 * (dis - s->rest_length) * a_b;
            s->m1->forces  += constra;
            s->m2->forces  += - constra;
        }

        for (auto &m : masses)
        {
            if (!m->pinned)
            {

                // m->position += v_next * delta_t;
                // m->velocity= v_next;
                Vector2D ac = ( gravity) / m->mass;
                Vector2D tmp_position = m->position;
                m->position +=  (1 - 0.00005)* (m->position - m->start_position) + ac * delta_t * delta_t + m->forces ;
                m->start_position = tmp_position;
                m->forces = Vector2D(0, 0);
                // TODO (Part 4): Add global Verlet damping
                
            }
        }
    }
}
