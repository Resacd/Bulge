#pragma once
#include "Node.h"
#include <fstream>
#include <string>
#include <mpi.h>
#include "Utility.h"
#include <filesystem>

#ifndef _WIN32
#include <sys/stat.h>
#else
#include <windows.h>
#endif

#include <memory> 

#include "Struct.h"
#include "HaloLayer.h"
#include "Subdomain.h"

#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPoints.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkDelaunay3D.h>
#include <vtkXMLUnstructuredGridWriter.h>

#include <vtkPolyData.h>
#include <vtkCellArray.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkFloatArray.h>
#include <vtkTriangleFilter.h>
#include <vtkPolyDataWriter.h>

#include <vtkHexahedron.h>
#include <vtkCellArray.h>
#include <vtkSmoothPolyDataFilter.h>
#include <vtkGeometryFilter.h>
#include <vtkXMLPolyDataWriter.h>



class OutputFile
{
private:

	std::mutex m_fileMutex;
	int m_outputIdx;
	std::mutex m_dataMutex;
	Utility Utility_;

public:

	OutputFile(int outputIdx);
	~OutputFile();


	void WriteVTK_LocalData(SubDomain &subdomain);
	void WriteWallShearStressVTK(SubDomain& subdomain, int iter);

};
