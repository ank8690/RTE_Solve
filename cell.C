#include"cell.H"

/*****************************Constructor Definitions*************************/
Cell::Cell()
{
	m_l_nbCells = 0;
	m_l_cellIndex = 0;
	m_l_nbRays = 0;
	m_l_vec_FaceIndicies.clear();
	m_d_vec_cell_ray_Inten.clear();
	/*if(m_d_vec_vec_cell_ray_Inten.size())
	{
		for(int i = 0; i < m_d_vec_vec_cell_ray_Inten.size(); ++i)
		{
			m_d_vec_vec_cell_ray_Inten[i].clear();
		}
		m_d_vec_vec_cell_ray_Inten.clear();
	}*/
}

/*****************************************************************************/
Cell::Cell(long l_nbCells)
{
	m_l_nbCells = l_nbCells;
	m_l_cellIndex = 0;
	m_l_nbRays = 0;
	m_l_vec_FaceIndicies.clear();
	m_d_vec_cell_ray_Inten.clear();
	/*if(m_d_vec_vec_cell_ray_Inten.size())
	{
		for(int i = 0; i < m_d_vec_vec_cell_ray_Inten.size(); ++i)
		{
			m_d_vec_vec_cell_ray_Inten[i].clear();
		}
		m_d_vec_vec_cell_ray_Inten.clear();
	}*/
}


/*****************************************************************************/
Cell::~Cell()
{
	m_l_nbCells = 0;
	m_l_cellIndex = 0;
	m_l_nbRays = 0;
	m_l_vec_FaceIndicies.clear();
	m_d_vec_cell_ray_Inten.clear();
	/*if(m_d_vec_vec_cell_ray_Inten.size())
	{
		for(int i = 0; i < m_d_vec_vec_cell_ray_Inten.size(); ++i)
		{
			m_d_vec_vec_cell_ray_Inten[i].clear();
		}
		m_d_vec_vec_cell_ray_Inten.clear();
	}*/
}
int main()
{
	return 1; 
}
