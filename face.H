/*--------------------------------------------------------------------------
Class
    Face
Description
    A Face is defined as a list of faces with extra functionality
SourceFiles
    cellI.H
---------------------------------------------------------------------------*/

#ifndef faceRAT_H
#define faceRAT_H
	
#include<vector>

//using namespace std;
class Face
{
	
			//long m_l_nbCells;
			long m_l_FaceIndex;
	public:
			
			long m_l_nbRays;
			//std::vector<long>m_l_vec_FaceIndicies;
			std::vector<double> m_d_vec_face_ray_Inten;	

			// constructors
			Face();
			//Face(long l_nbCells);

			// destructor
			~Face();
			
			
};
#endif
