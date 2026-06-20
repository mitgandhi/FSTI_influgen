#include "fem.h"
#include "string_functions.h"
#include <map>

using namespace std;

void fem::setupgap(std::string gap_set_name, im_direction_t direction)
{
	//set the gap set name
	if(nodesets.count(gap_set_name) == 0)
	{
		cout << "ERROR: The node set named \"" << gap_set_name << "\" is required but not defined in the mesh!" << endl;
		exit(1);
	}
	gap_nodeset = &nodesets[gap_set_name];

	//set the im normal direction
	im_direction = direction;
	
	//create the gap face set
	nodeset2faceset(*gap_nodeset, gap_face_vector);
	cout << "  Gap face set created with " << gap_face_vector.size() << " faces." << endl;

	//make the gap_face_set for from the gap_face_vector to ensure other set loads don't also load the gap
	for(auto f=gap_face_vector.begin(); f!=gap_face_vector.end(); ++f)
	{
		gap_face_set.insert(*f);
	}
	
	//create the gap node norms (actually an average of attached face norms)
	gap_node_norms.resize(gap_nodeset->size());
	std::vector<double> gap_node_norms_cnt(gap_nodeset->size(), 0);

	// build a map of nodes -> faces
	std::map<const node *, std::vector<const face *> > gapnode2gapfaces;
	for(auto f=gap_face_set.begin(); f!=gap_face_set.end(); ++f)
	{
		for(size_t n=0; n<(*f)->nodes.size(); ++n)
		{
			gapnode2gapfaces[(*f)->nodes[n]].push_back(*f);
		}
	}

	//now loop over the gap_nodeset and calc the norm for the faces attached to each node
	for(size_t n=0; n<gap_nodeset->size(); ++n)
	{
		const node* N = &nodes[(*gap_nodeset)[n]];
		const std::vector<const face *>& Faces = gapnode2gapfaces[N];
		//loop over the node faces
		for(size_t f=0; f<Faces.size(); ++f)
		{
			gap_node_norms[n] = gap_node_norms[n]+Faces[f]->normal();
			gap_node_norms_cnt[n] += 1;
		}
	}

	//now do the division by gap_node_norms_cnt
	for(size_t n=0; n<gap_nodeset->size(); ++n)
	{
		gap_node_norms[n] /= gap_node_norms_cnt[n];
	}
	
}
void fem::nodeset2faceset(const nodeset& ns, std::vector<const face*>& fs)
{
	//build a std::set of node* from the nodeset for quick searching
	std::set<const node*> ns_set;

	for(auto n=ns.begin(); n!=ns.end(); ++n)
	{
		ns_set.insert(&nodes[*n]);
	}

	//used below, but ensuring for optimization that this comparison is only done once
	const bool is_gap_nodeset = ns == *gap_nodeset;

	//now loop over all the faces
	for(auto f=externalelementfaces.begin(); f!=externalelementfaces.end(); ++f)
	{
		bool set_face = true;
		for(size_t n=0; n<f->first.nodes.size(); ++n)
		{
			if(ns_set.count(f->first.nodes[n]) == 0)	//the face node isn't in the nodeset
			{
				set_face = false;
			}
		}
		if(set_face)
		{
			//don't insert a face if it is part of the gap face set
			//unless this nodeset is the gap node set
			if(gap_face_set.count(&(f->first)) == 0 || is_gap_nodeset)
			{
				fs.push_back(&(f->first));
			}
		}
	}
}
void fem::load_gap_face(size_t faceid, double pressure, loadmap_t& loadmap)
{
	const face& f = *gap_face_vector[faceid];
	double area = f.area();
	point nodeload = f.normal();
	nodeload *= pressure*area/double(f.nodes.size());
	for(size_t n=0; n<f.nodes.size(); ++n)
	{
		loadmap[f.nodes[n]->nid] += nodeload;
	}
}
void fem::load_set(std::string setname, double pressure, loadmap_t& loadmap)
{
	setname = string_functions::str2lower(setname);

	if(nodesets.count(setname) != 1)
	{
		cout << "ERROR: Unable to load the set \"" << setname << "\" becuase it does not exist in the mesh!" << endl;
		exit(1);
	}
	//first convert the nodeset to a faceset
	std::vector<const face*> faceset;
	nodeset2faceset(nodesets[setname], faceset);

	//now loop over all the set faces
	for(auto f_it=faceset.begin(); f_it!=faceset.end(); ++f_it)
	{	
		const face& f = **f_it;
		double area = f.area();
		point nodeload = f.normal();
		nodeload *= pressure*area/double(f.nodes.size());
		for(size_t n=0; n<f.nodes.size(); ++n)
		{
			loadmap[f.nodes[n]->nid] += nodeload;
		}
	}
}
std::vector<double> fem::inverse_3x3(const vector<double> A)
{
	vector<double> B(9, 0);

	//det(A)
	double det = A[0]*( A[4]*A[8]-A[7]*A[5] )
                 - A[1]*( A[3]*A[8]-A[5]*A[6] )
                 + A[2]*( A[3]*A[7]-A[4]*A[6] );
	
	//Inverse of A
	double invdet = 1.0/det + 1.0e-12;
	B[0] =  ( A[4]*A[8] - A[7]*A[5] ) * invdet;
	B[3] = -( A[1]*A[8] - A[2]*A[7] ) * invdet;
	B[6] =  ( A[1]*A[5] - A[2]*A[4] ) * invdet;
	B[1] = -( A[3]*A[8] - A[5]*A[6] ) * invdet;
	B[4] =  ( A[0]*A[8] - A[2]*A[6] ) * invdet;
	B[7] = -( A[0]*A[5] - A[3]*A[2] ) * invdet;
	B[2] =  ( A[3]*A[7] - A[6]*A[4] ) * invdet;
	B[5] = -( A[0]*A[7] - A[6]*A[1] ) * invdet;
	B[8] =  ( A[0]*A[4] - A[3]*A[1] ) * invdet;

	return B;
};
void fem::loadmap2b(const loadmap_t& loadmap, std::vector<double> &b)
{
	//initialize the load vector
	b.clear();
	b.resize(udof, 0);

	loadmap2b(loadmap, b.begin());
}
void fem::loadmap2b(const loadmap_t& loadmap, std::vector<double>::iterator b)
{
	//this assumes that the vector underlying the b iterator is AT LEAST b+udof long
	
	//set the loads from the map into a b vector
	for(auto n=loadmap.begin(); n!=loadmap.end(); ++n)
	{
		for(int dof=0; dof<3; ++dof)
		{
			int gdof = nodes[n->first].DOF[dof];
			if(gdof >= 0)	//node is unconstrained
			{
				*(b+gdof) += n->second.coord[dof];
			}
		}
	}
}
size_t fem::get_gap_facecnt()
{
	return gap_face_vector.size();
}
int fem::get_udof()
{
	return udof;
}
size_t fem::get_gap_nodecnt()
{
	return gap_nodeset->size();
}
