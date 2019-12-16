#include "MatLabVariable.h"
#include <array>
//#include <boost/locale.hpp>
//#include <boost/format.hpp>

//#include "../../../Exception.h"
//#include "../../../StringUtilities.h"

#define var const auto

using namespace Eigen;
namespace Jde::IO::MatLab
{
#pragma region Constructors
	MatLabVariable::MatLabVariable( string_view name, const Eigen::MatrixXd& matrix ):
		_name{name},
		_matrixType(Matrix),
		_matrix(&matrix)
	{}

	MatLabVariable::MatLabVariable( string_view name, const Eigen::VectorXd& vector ):
		_name{name},
		_matrixType(Vector),
		_vector(&vector)
	{}
	MatLabVariable::MatLabVariable( string_view name, const std::shared_ptr<Eigen::SparseMatrix<double>> matrix ):
		_name{name},
		_matrixType(Sparse),
		_sparse(matrix)
	{}

	MatLabVariable::MatLabVariable( string_view name, matio_classes classType, matio_types data_type, int rank, size_t *dims, void *data, int opt ):
		_pVariable( Mat_VarCreate(string(name).c_str(), classType, data_type, rank, dims, data, opt) ),
		_name( name ),
		_matrixType(Data)
	{}

	MatLabVariable::MatLabVariable( matvar_t* variable ):
		_pVariable(variable),
		_name( variable->name )
	{}

	MatLabVariable::MatLabVariable( const MatLabVariable& copyFrom ):
		_pVariable( copyFrom._pVariable),
		_name(copyFrom._name),
		_matrixType(copyFrom._matrixType)
	{
		if( _matrixType==Matrix )
			_matrix = copyFrom._matrix;
		else if( _matrixType==Vector )
			_vector= copyFrom._vector;
		else if( _matrixType==Sparse )
			_sparse = copyFrom._sparse;
	}

	MatLabVariable::~MatLabVariable()
	{ 
		//Logging::ProcessInfo::LogMemoryUsage( (boost::wformat("Before Variable '%1%' Free") % _name).str() );
		Mat_VarFree(_pVariable);
		_pVariable = nullptr;
		//Logging::ProcessInfo::LogMemoryUsage( (boost::wformat("After Variable '%1%' Free") % _name).str() );
	}
#pragma endregion
/*
	sp<Eigen::MatrixXd> MatLabVariable::ToMatrix()const
	{
		//MatLab::MatLabVariable variable = file[variableName];
		if( _pVariable==NULL )
			throw Exception( "Variable not initialized" );
		
		size_t stride = Mat_SizeOf( _pVariable->data_type );
		
		if( _pVariable->class_type!=MAT_C_DOUBLE )
			throw Exception( "class_type!=MAT_C_DOUBLE" );
		
		const size_t rowCount =  _pVariable->dims[0]; 
		const size_t columnCount =  _pVariable->dims[1];
		auto matrix = make_shared<Eigen::MatrixXd>( rowCount,columnCount );
		for( size_t i = 0; i < _pVariable->dims[0]; ++i ) 
		{
			for ( size_t j = 0; j < _pVariable->dims[1]; ++j ) 
			{
				size_t index = _pVariable->dims[0]*j+i;
				const char* data = static_cast<const char*>(_pVariable->data);
				double* value = (double*)( data+index*stride );
				(*matrix)(i,j) = *value;
			}
		}
		return matrix;
	}
	sp<Eigen::VectorXd> MatLabVariable::ToVector()const
	{
		if( _pVariable==NULL )
			throw Exception( "Variable not initialized" );

		size_t stride = Mat_SizeOf( _pVariable->data_type );

		if( _pVariable->class_type!=MAT_C_DOUBLE )
			throw Exception( "class_type!=MAT_C_DOUBLE" );

		const size_t rowCount =  _pVariable->dims[0]; 
		//const size_t columnCount =  _pVariable->dims[1];
		auto vector = make_shared<Eigen::VectorXd>(rowCount);	
		for( size_t index = 0; index < rowCount; ++index ) 
		{
			const char* data = static_cast<const char*>(_pVariable->data);
			double* value = (double*)( data+index*stride );
			(*vector)(index) = *value;
		}
		return vector;
	}
*/
	double MatLabVariable::operator()(size_t rowIndex, size_t columnIndex) const
	{
		const char* data = static_cast<const char*>(_pVariable->data);
		auto index = _pVariable->dims[0]*columnIndex+rowIndex;
		auto stride = Mat_SizeOf( _pVariable->data_type );
		double* value = (double*)( data+index*stride );
		return *value;
	}
	void MatLabVariable::Write( mat_t* pFile )const
	{
		if( _matrixType!=Data )
		{
			if( _matrixType==Sparse )
			{
				Write( pFile, _name, *_sparse );
				return;
			}
			else
			{
				const size_t rowCount = _matrixType==Matrix ? _matrix->rows() : _vector->rows();
				const size_t columnCount = _matrixType==Matrix ? _matrix->cols() : 1;
				size_t* dims = static_cast<size_t*>( malloc( 2*sizeof(size_t)) ); 
				dims[0] = rowCount;
				dims[1] = columnCount;
				double* values = static_cast<double*>( malloc(rowCount*columnCount*sizeof(double)) ); 
				for( uint rowIndex=0; rowIndex<rowCount; ++rowIndex )
				{
					for( uint columnIndex=0; columnIndex<columnCount; ++columnIndex )
						values[rowIndex*columnCount+columnIndex] = _matrixType==Matrix ? (*_matrix)(rowIndex, columnIndex) : (*_vector)(rowIndex);
				}
				_pVariable = Mat_VarCreate( _name.c_str(), MAT_C_DOUBLE,MAT_T_DOUBLE,2,dims,values,0);
			}
		}
		Mat_VarWrite( pFile, _pVariable, MAT_COMPRESSION_ZLIB/*MAT_COMPRESSION_NONE*/ );
	}
#pragma region GetMatioType
	matio_types MatLabVariable::GetMatioType( const std::type_info& typeInfo )
	{
		matio_types dataType = matio_types::MAT_T_DOUBLE;
		if( typeInfo==typeid(float) )
			dataType = matio_types::MAT_T_SINGLE;
		else if( typeInfo==typeid(int) )
			dataType = matio_types::MAT_T_INT32;
		else if( typeInfo==typeid(size_t) )
			dataType = matio_types::MAT_T_UINT64;
		else if( typeInfo!=typeid(double) )
			THROW( Exception( fmt::format( "Datatype '{}' has not been implemented.", typeInfo.name())) );
		return dataType;
	}
	matio_classes MatLabVariable::GetMatioClass( const std::type_info& typeInfo )noexcept(false)
	{
		matio_classes classType = matio_classes::MAT_C_DOUBLE;
		if( typeInfo==typeid(float) )
			classType = matio_classes::MAT_C_SINGLE;
		else if( typeInfo==typeid(int) )
			classType = matio_classes::MAT_C_INT32;
		else if( typeInfo==typeid(size_t) )
			classType = matio_classes::MAT_C_UINT64;
		else if( typeInfo!=typeid(double) )
			THROW( Exception(fmt::format("class '{}' has not been implemented.", typeInfo.name())) );
		return classType;
	}
#pragma endregion
	const mat_sparse_t* MatLabVariable::GetSparsePointer()const noexcept(false)
	{
		if( _pVariable->class_type!=MAT_C_SPARSE )
			THROW( Exception( "class_type!=MAT_C_SPARSE") );
		return static_cast<mat_sparse_t*>(_pVariable->data);
	}


	const std::vector<std::pair<size_t,size_t>> MatLabVariable::ColumnCounts()const noexcept(false)
	{
		std::vector<std::pair<size_t,size_t>> columnCounts( cols() );
		if( _pVariable->class_type==MAT_C_SPARSE )
		{
			const mat_sparse_t* sparse =  GetSparsePointer();
			size_t startIndex=0;
			for( mat_uint32_t jcIndex=1; jcIndex<sparse->njc; ++jcIndex )
			{
				const auto currentIndex = sparse->jc[jcIndex];
				std::pair<size_t,size_t> value = make_pair( currentIndex-startIndex,startIndex );
				columnCounts[jcIndex-1] = value;
				startIndex = value.first+value.second;
			}
		}
		else
		{
			for( uint columnIndex=0; columnIndex<cols(); ++columnIndex )
				columnCounts[columnIndex] = make_pair( rows(), rows()*columnIndex );
		}
		return columnCounts;
	}

	std::vector<string> WToColumnNames( matvar_t& variable )
	{
		const size_t rowCount = variable.dims[0];
		const size_t columnCount = variable.dims[1];
		std::vector<string> columnNames( rowCount );
		const char* pData2 = static_cast<char*>( variable.data );
		//const uint16_t* pData = static_cast<const uint16_t*>( variable.data );
		for( uint columnIndex=0, dataIndex=0; columnIndex<columnCount; ++columnIndex )
		{
		//	DBG( "columnIndex={}", columnIndex );
			for( uint rowIndex=0; rowIndex<rowCount; ++rowIndex, dataIndex+=2 )
			{
				auto& value = columnNames[rowIndex];
				var ch = pData2[dataIndex];
				value += static_cast<char>( ch );
			//	if( columnIndex<3 )
			//		DBG( "{}={}, {}", rowIndex, value, dataIndex );
			}
		}
		for( auto& columnName : columnNames )
			StringUtilities::RTrim( columnName );

		return columnNames;
	}
	std::vector<string> MatLabVariable::ToColumnNames()
	{
		if( _pVariable->data_type==MAT_T_UINT16 )
			return WToColumnNames( *_pVariable );
		const size_t rowCount = _pVariable->dims[0];
		const size_t columnCount = _pVariable->dims[1];
		std::vector<string> columnNames( rowCount );
		const char* pData = static_cast<const char*>( _pVariable->data );
		for( uint columnIndex=0, dataIndex=0; columnIndex<columnCount; ++columnIndex )
		{
			for( uint rowIndex=0; rowIndex<rowCount; ++rowIndex )
			{
				auto& value = columnNames[rowIndex];
				value += pData[dataIndex++];
			}
		}
		for( auto& columnName : columnNames )
			StringUtilities::RTrim( columnName );

		return columnNames;
	}
}