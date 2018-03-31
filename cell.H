/*--------------------------------------------------------------------------
Class
    Cell
Description
    A cell is defined as a list of faces with extra functionality
SourceFiles
    cellI.H
---------------------------------------------------------------------------*/

#ifndef cellRAT_H
#define cellRAT_H

#include<vector>
class Cell
{
	
			long m_l_nbCells;
			long m_l_cellIndex;
	public:
			long m_l_FaceIndex_East;
			long m_l_FaceIndex_West;
			long m_l_FaceIndex_North;
			long m_l_FaceIndex_South;
			long m_l_nbRays;
			float m_f_Vol;
			std::vector<long>m_l_vec_FaceIndicies;
			std::vector<double> m_d_vec_cell_ray_Inten;	

			// constructors
			Cell();
			Cell(long l_nbCells);

			// destructor
			~Cell();
			
			
};
#endif
