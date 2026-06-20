#include "fem.h"
#include "dmatrix.h"
#include "element_C3D4.h"
using namespace std;

void fem::number_DOF()
{
	//we will first reset all the node DOF's
	//then mark the constrained DOF's
	//then renumber all but the constrainted DOF's
	for(size_t n=0; n<nodes.size(); ++n)
	{
		const size_t ndof = nodes[n].DOF.size();
		for(size_t dof=0; dof<ndof; ++dof)
		{
			nodes[n].DOF[dof] = 0;
		}
	}

	if(!options.inrel)	//ignore constraints if IR is enabled
	{
		for(auto nset=nodesets.begin(); nset!=nodesets.end(); ++nset)
		{
			//only bother looping over nodesets where a constraint is set
			if(nset->second.const_x || nset->second.const_y || nset->second.const_z)	
			{
				for(auto n=nset->second.begin(); n!=nset->second.end(); ++n)
				{
					if(nset->second.const_x)
					{
						nodes[*n].DOF[0] = -1;
					}
					if(nset->second.const_y)
					{
						nodes[*n].DOF[1] = -1;
					}
					if(nset->second.const_z)
					{
						nodes[*n].DOF[2] = -1;
					}
				}
			}
		}
	}

	udof = 0;
	for(size_t n=0; n<nodes.size(); ++n)
	{
		const size_t ndof = nodes[n].DOF.size();
		for(size_t dof=0; dof<ndof; ++dof)
		{
			if(nodes[n].DOF[dof] != -1)
			{
				nodes[n].DOF[dof] = udof;
				++udof;
			}
		}
	}

	//add 3 DOF in the case of IR for the lambda constraints
	if(options.inrel)
	{
		udof += 6;
	}
};
std::vector<const node*> fem::nid2node(const std::vector<size_t> &nids)
{
	const size_t ncnt = nids.size();
	vector<const node*> n(ncnt);
	
	for(size_t i=0; i<ncnt; ++i)
	{
		n[i] = &nodes[nids[i]];
	}
	return n;
}

void fem::stiffness_matrix(const bool only_triu)
{
	//first number the DOF's
	number_DOF();
	cout << "Number of unconstrainted DOF: " << udof << endl;

	//Resize the K matrix
	K_gmm.clear_mat();
	K_gmm.resize(udof, udof);
	nodes_mass.clear();
	nodes_mass.resize(nodecnt, 0);

	cout << "Building stiffness matrix:" << endl;
	int eprocessed = 0;	//used for progress status
	double lastpdisp = 0;	//used for cout'ing progress

	//the average (abs) stiffness K value will be used in case of IR for the constraint
	double avgK = 0;
	int avgK_cnt = 0;

	//loop over all the elements, get the local stiffness matrix
	//assemble it into the global matrix
	//also store the nodal mass, to be used in the case of IR

	//we actually have to loop over the elementsets b/c that is where materials are defined
	for(auto eset_ptr=elementsets.begin(); eset_ptr!=elementsets.end(); ++eset_ptr)
	{
		elementset* eset = &eset_ptr->second;	//a pointer to the element set

		//now for every element in the elementset
		for(int eset_e=0; eset_e<eset->size(); ++eset_e)
		{
			size_t global_e = (*eset)[eset_e];	//the element id in the global element vector

			dmatrix<double> Kele;	//the local element stiffness matrix
			vector<double> ele_node_mass;	//the local element node mass vector
			vector<face> ele_faces;	//vector of element faces
			vector<int> ldof2gdof;	//a vector to map the local element dof 2 the global dof

			//here is where more element types could easily be added
			switch(elements[global_e].type)
			{
			case element::C3D4:
					C3D4 ele(nid2node(elements[global_e].nodes), eset->material);
					ele.calc(Kele, ele_node_mass);
					ldof2gdof = ele.ldof2gdof();
					ele_faces = ele.getfaces();
				break;
			}
			
			//add the local element stiffness matrix to the global matrix
			for(int i=0; i<Kele.m; i++)
			{
				for(int j=0; j<Kele.n; j++)
				{
					if(ldof2gdof[i] >= 0 && ldof2gdof[j] >= 0)
					{
						if(ldof2gdof[j] >= ldof2gdof[i] || !only_triu)	//used to only store the upper half
						{
							K_gmm(ldof2gdof[i], ldof2gdof[j]) += Kele(i,j);
							avgK += fabs(Kele(i,j));
							++avgK_cnt;
						}
					} else {
						//node is constrained

						//if non-zero essential boundary conditions are specified, we would have to do more here,
						//but we don't since only 0 val is currently supported.
					}
				}
			}

			//add the element faces, if the face is already found, it is not an external face
			//so remove it
			for(size_t i=0; i<ele_faces.size(); ++i)
			{
				auto f = externalelementfaces.find(ele_faces[i]);
				if(f == externalelementfaces.end())
				{
					externalelementfaces[ele_faces[i]] = global_e;
				} else {
					externalelementfaces.erase(f);
				}
			}

			//add the element node mass to the global node mass vector
			for(size_t i=0; i<ele_node_mass.size(); ++i)
			{
				nodes_mass[elements[global_e].nodes[i]] += ele_node_mass[i];
			}
				
			//update the status counter
			eprocessed++;

			//used to display progress update
			double p = double(eprocessed)/double(elementcnt);
			if(p  > lastpdisp + 0.05)
			{
				cout << "  Progress: " << int(p*100.0) << "%" << endl;
				lastpdisp = p;
			}
		}
	}


	if(options.inrel == 1)
	{
		//loop over all the nodes to build the body mass, cog, and avg nodal mass
		IR.body_mass = 0;
		IR.cog = point(0,0,0);
				
		for(size_t n=0; n<nodecnt; ++n)
		{
			IR.body_mass += nodes_mass[n];
			IR.cog += nodes[n].coord*nodes_mass[n];
		}
		IR.cog /= IR.body_mass;
		double avgNodeMass = IR.body_mass/double(nodecnt);

		//Scale factor = [abs mean stiffness matrix] / [mean node mass]
		//This significantly improves the condition number of the stiffness matrix
		double scale = (avgK/double(avgK_cnt))/avgNodeMass;
		IR.scale = scale;

		//ax, ay, az const
		for(int n=0; n<nodecnt; n++)
		{
			for(int d=0; d<3; d++)
			{
				const int Gdof = nodes[n].DOF[d];

				K_gmm(Gdof, udof - 3 + d) = scale*nodes_mass[n];
				if(!only_triu)	//used to store the lower half when necessary
				{
					K_gmm(udof - 3 + d, Gdof) = K_gmm(Gdof, udof - 3 + d);
				}
			}
		}

		//alphax, alphay, alphaz
		for(int n=0; n<nodecnt; n++)
		{
			double rx = nodes[n].coord.x()-IR.cog.x();
			double ry = nodes[n].coord.y()-IR.cog.y();
			double rz = nodes[n].coord.z()-IR.cog.z();
			double mn = scale*nodes_mass[n]/(rx*rx+ry*ry+rz*rz);

			//Fx = ((rx*rx+rz*rz)*rz+ry*ry*rz)*ay+(-ry*rz*rz-(rx*rx+ry*ry)*ry)*az
			{
				const int Gdof = nodes[n].DOF[0];
				K_gmm(Gdof, udof-5) = mn*((rx*rx+rz*rz)*rz+ry*ry*rz);
				K_gmm(Gdof, udof-4) = mn*(-ry*rz*rz-(rx*rx+ry*ry)*ry);
			}

			//Fy = (-(ry*ry+rz*rz)*rz-rx*rx*rz)*ax+(rx*rz*rz+(rx*rx+ry*ry)*rx)*az
			{
				const int Gdof = nodes[n].DOF[1];
				K_gmm(Gdof, udof-6) = mn*(-(ry*ry+rz*rz)*rz-rx*rx*rz);
				K_gmm(Gdof, udof-4) = mn*(rx*rz*rz+(rx*rx+ry*ry)*rx);
			}

			//Fz = ((ry*ry+rz*rz)*ry+rx*rx*ry)*ax+(-rx*ry*ry-(rx*rx+rz*rz)*rx)*ay
			{
				const int Gdof = nodes[n].DOF[2];
				K_gmm(Gdof, udof-6) = mn*((ry*ry+rz*rz)*ry+rx*rx*ry);
				K_gmm(Gdof, udof-5) = mn*(-rx*ry*ry-(rx*rx+rz*rz)*rx);
			}

			if(!only_triu)	//used to store the lower half when necessary
			{
				for(int d=0; d<3; d++)
				{
					const int Gdof = nodes[n].DOF[d];
					K_gmm(udof-6, Gdof) = K_gmm(Gdof, udof-6);
					K_gmm(udof-5, Gdof) = K_gmm(Gdof, udof-5);
					K_gmm(udof-4, Gdof) = K_gmm(Gdof, udof-4);
				}		
			}	
		}
	}

	//cout finished progress
	cout << "  Finished!" << endl;
};

