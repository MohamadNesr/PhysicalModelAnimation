/*
 * CalculsMSS.cpp :
 * Copyright (C) 2016 Florence Zara, LIRIS
 *               florence.zara@liris.univ-lyon1.fr
 *               http://liris.cnrs.fr/florence.zara/
 *
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 2.1 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

/** \file CalculsMSS.cpp
Programme calculant pour chaque particule i d un MSS son etat au pas de temps suivant 
 (methode d 'Euler semi-implicite) : principales fonctions de calculs.
\brief Fonctions de calculs de la methode semi-implicite sur un systeme masses-ressorts.
*/ 

#include <stdio.h>
#include <math.h>
#include <vector>
#include <iostream>

#include "vec.h"
#include "ObjetSimule.h"
#include "ObjetSimuleMSS.h"
#include "Viewer.h"

using namespace std;





/**
* Calcul des forces appliquees sur les particules du systeme masses-ressorts.
 */
void ObjetSimuleMSS::CalculForceSpring()
{
	/// f = somme_i (ki * (l(i,j)-l_0(i,j)) * uij ) + (nuij * (vi - vj) * uij) + (m*g) + force_ext
	
	/// Rq : Les forces dues a la gravite et au vent sont ajoutees lors du calcul de l acceleration
	    for (auto i = 0; i < _SystemeMasseRessort->GetNbParticule(); i++){
			
			Particule *currentParticule = _SystemeMasseRessort->GetPartList()[i];
			Vector somme_i = Vector(0.0, 0.0, 0.0);
			for(auto j=0; j < currentParticule->GetNbVoisins(); j++){

				Ressort *r = currentParticule->GetRessortList()[j];
				Particule *particuleA = r->GetParticuleA();
				Particule *particuleB = r->GetParticuleB();
				//si la particule courante est B on inverse
				if (particuleB == _SystemeMasseRessort->GetParticule(i)){
					std::swap(particuleA,particuleB);
				}

				// ki
				float ki = r->GetSpring()->_Raideur;
				// l(i,j)
				float l_ij = distance(Point(P[particuleA->_Id]), Point(P[particuleB->_Id]));
				// l_0(i,j))
				float l_0_ij = r->GetSpring()->_L0;
				// uij
				Vector uij = normalize(P[particuleB->_Id] - P[particuleA->_Id]);
				// (ki * (l(i,j)-l_0(i,j)) * uij )
				Vector kllu = ki * (l_ij - l_0_ij) * uij;

				Vector AB = particuleB->GetPosition() - particuleA->GetPosition();
				// vecteur
				Vector v = AB/l_ij;
				// amortissement : F = -coef_amortissemnet.vecteur
				float f_amortissemnet = dot(-r->GetAmortissement(), v);

				// on applique l'amortissement a nos ressorts
				Vector dynamique = AB/length(AB)*f_amortissemnet;
				
				kllu = kllu + dynamique;
				somme_i = somme_i + kllu;

				// dechirure
				if(l_ij>1){
					currentParticule->GetRessortList().erase(currentParticule->GetRessortList().begin()+j);
				}
				
			}
			Force[currentParticule->_Id] = somme_i;
		}		
}//void
 //pt fixe 0 70

/**
 * Gestion des collisions avec le sol.
 */
void ObjetSimuleMSS::Collision()
{
    /// Arret de la vitesse quand touche le plan
    for(auto i=0; i<P.size(); i++){
        if(P[i].y < -9.99){
			P[i].y = -9.99;
        	V[i] = 0;
       	}
    }	
    
}// void

/**
 * Gestion des collisions avec la sphere.
 */
void ObjetSimuleMSS::CollisionSphere(Point p, double rayon) {
    for (auto i = 0; i < P.size(); i++) {
        if( distance(Point(P[i]), p) < rayon) {   
 			Vector cp(p,Point(P[i]));
            Vector n = normalize(cp);
            P[i] = Vector(n * rayon + p);
			V[i]= 0;
        }
    }
}// void