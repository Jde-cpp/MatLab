#include "MatLabFile.h"
#include <array>
//#include <boost/format.hpp>
#include <mutex>
#include <tuple>

//#include "../../../Exception.h"
#include "MatLabVariable.h"
//#include "../../../Stopwatch.h"
#define var const auto

namespace Jde::IO::MatLab
{
	std::mutex* MatLabFile::_lock = new std::mutex();
	std::list<std::tuple<int,std::string>>* MatLabFile::_logs = nullptr;
#pragma region Constructors
	MatLabFile::MatLabFile()
	{}
	MatLabFile::MatLabFile( const fs::path& matlabFileName )noexcept(false):
		_fileName(matlabFileName)
	{
		//boost::format format = ;
		//auto error = boost::str( boost::format{"(%1%) Could not open file:  '%2%'"} % 5 % "foo" );
		//auto error = boost::str( boost::format{"%2%/%1%/%3%"} % 12 % 5 % 2014 );
		if( _fileName.string().size()==0 )
			THROW( Exception("matlabFileName is empty"sv) );
		_lock->lock();
		if( _logs==nullptr )
		{
			Mat_LogInitFunc( "Jde.Math.Matlab.MatLabFile", MatLabFile::Logfunc );
			_logs = new std::list<std::tuple<int,std::string>>();
		}
		else
			_logs->clear();
		_file = Mat_Open( _fileName.string().c_str(), MAT_ACC_RDONLY );
		if( _file==nullptr )
		{
			var error =  _logs->size()>0  ? std::get<1>( *_logs->begin() ) : std::to_string(errno);
			THROW( Exception( "({}) Could not open file:  '{}'"sv, error, _fileName) );
			//throw std::exception( error2.c_str() );
		}
		_lock->unlock();
	}

	MatLabFile::~MatLabFile()
	{
		if( _file )
			Mat_Close( _file );
	}
#pragma endregion


#pragma region operator[]
/*	MatLabVariable MatLabFile::operator[]( const wchar_t* const variableName )
	{
		return (*this)[boost::locale::conv::utf_to_utf<char>(variableName).c_str()];
	}*/
	MatLabVariable MatLabFile::operator[]( string_view variableName )
	{
		//Stopwatch sw( StopwatchTypes::ReadFile, variableName, true );
		matvar_t* pVariable = Mat_VarRead( _file, string(variableName).c_str() );
		if( !pVariable )
			THROW( Exception("Could not find variable '{}' in file '{}'"sv, variableName, _fileName) );
		return MatLabVariable( pVariable );
	}
#pragma endregion
	void MatLabFile::Logfunc( int log_level, char *message )
	{
		auto value = std::make_tuple( log_level, std::string(message) );
		_logs->push_back( value );
	}
	void MatLabFile::Combine( string_view path, string_view variableName, MatLabVariable v1, MatLabVariable v2 )
	{
		Stopwatch sw( fmt::format("StopwatchTypes::WriteFile - {}", variableName), true );

		matvar_t* pVariable1 = v1._pVariable;
		//const size_t stride = Mat_SizeOf( pVariable1->data_type );
		if( pVariable1->class_type!=MAT_C_SPARSE )
			THROW( Exception("class_type!=MAT_C_SPARSE"sv) );
		matvar_t* pVariable2 = v2._pVariable;

		mat_sparse_t* pSparse1 =  static_cast<mat_sparse_t*>(pVariable1->data);
		mat_sparse_t* pSparse2 =  static_cast<mat_sparse_t*>(pVariable2->data);

		mat_sparse_t* pCombined = static_cast<mat_sparse_t*>( calloc(1, sizeof(mat_sparse_t)) );

		pCombined->nir = pCombined->ndata = pCombined->nzmax = pSparse1->nzmax+pSparse2->nzmax;

		mat_uint32_t* pRowIndexes = static_cast<mat_uint32_t*>( calloc( pCombined->nir,sizeof(mat_uint32_t)) );
		DBG( "pRowIndexes:  {:x}.  size:  {}"sv, (size_t)pRowIndexes, pCombined->nir*sizeof(int) );
		if( pRowIndexes==nullptr )
			THROW( Exception("out of memory"sv) );
		DBG( "copy to pRowIndexes:  {}"sv, pSparse1->nir*sizeof(int) );
		std::copy( pSparse1->ir, pSparse1->ir+pSparse1->nir, pRowIndexes );
		DBG( "pRowIndexes:  {:x}, size:  {}"sv, (size_t)(pRowIndexes+pSparse1->nir*sizeof(int)), pSparse2->nir*sizeof(int) );
		std::copy( pSparse2->ir, pSparse2->ir+pSparse2->nir, pRowIndexes+pSparse1->nir );
		//memcpy( static_cast<void*>(pRowIndexes+pSparse1->nir*sizeof(int)), static_cast<void*>(pSparse2->ir), pSparse2->nir*sizeof(int) );
		pCombined->ir = pRowIndexes;
		DBG( "first:  {}.  vs:  {}  last:  {}  vs:  {}"sv, pRowIndexes[0], pSparse1->ir[0], pRowIndexes[pCombined->nir-1], pSparse2->ir[pSparse2->nir-1] );
		//	last:   (size_t)(pRowIndexes+pSparse1->nir*sizeof(int)) << std::dec << "size:  " << pSparse2->nir*sizeof(int) << endl;

		//const auto columnCount1 = pSparse1->njc-1;
		const auto columnCount = pSparse1->njc-1 + pSparse2->njc-1;
		pCombined->njc = columnCount+1;
		mat_uint32_t* jc = static_cast<mat_uint32_t*>( calloc(columnCount+1, sizeof(mat_uint32_t)) );
		for( uint jcIndex=1; jcIndex<pSparse1->njc; ++jcIndex )
		{
			jc[jcIndex] = pSparse1->jc[jcIndex];
			DBG( "jcIndex:  {} ({});  difference:  {}"sv, jcIndex, jc[jcIndex], jc[jcIndex]-jc[jcIndex-1] );
		}

		const auto startIndexCount=jc[pSparse1->njc-1];
		for( uint jcIndex=1; jcIndex<pSparse2->njc; ++jcIndex )
		{
			jc[pSparse1->njc+jcIndex-1] = startIndexCount+pSparse2->jc[jcIndex];
			DBG( "jcIndex:  {} ({});  difference:  {}"sv, pSparse1->njc+jcIndex-1, jc[pSparse1->njc+jcIndex-1], jc[pSparse1->njc+jcIndex-1]-jc[pSparse1->njc+jcIndex-2] );
		}
		pCombined->jc = jc;

		float* pValues = static_cast<float*>( calloc(pCombined->ndata,sizeof(float)) );
		float* pData1 = static_cast<float*>(pSparse1->data);
		std::copy( pData1, pData1+pSparse1->ndata, pValues );
		float* pData2 = static_cast<float*>(pSparse2->data);
		std::copy( pData2, pData2+pSparse2->ndata, pValues+pSparse1->ndata );
		//memcpy( static_cast<void*>(values+pSparse1->ndata*sizeof(float)), static_cast<void*>(pSparse2->data), pSparse2->ndata*sizeof(float) );
		DBG( "first:  {}.  vs:   {}.  last:  {}.  vs:  {}"sv, pValues[0], pData1[0], pValues[pCombined->ndata-1], pData2[pSparse2->ndata-1] );


		pCombined->data = pValues;

		size_t* dims = static_cast<size_t*>( malloc( 2*sizeof(size_t)) );
		const int rowCount = int(pVariable1->dims[0]);
		dims[0] = rowCount;
		dims[1] = columnCount;
		auto pVariable = Mat_VarCreate( string(variableName).c_str(), MAT_C_SPARSE,MatLabVariable::GetMatioType( typeid(float) ),2,dims,static_cast<void*>(pCombined), MAT_F_DONT_COPY_DATA );
		_file = Mat_CreateVer( string(path).c_str(), NULL, MAT_FT_MAT73 );
		Mat_VarWrite( _file, pVariable, MAT_COMPRESSION_ZLIB );
		{
			mat_sparse_t *sparse2 = (mat_sparse_t*)pVariable->data;
			if ( sparse2->ir != NULL )
				free(sparse2->ir);
			if ( sparse2->jc != NULL )
				free(sparse2->jc);
			if ( pVariable->isComplex && NULL != sparse2->data ) {
				mat_complex_split_t *complex_data = (mat_complex_split_t*)sparse2->data;
				free(complex_data->Re);
				free(complex_data->Im);
				free(complex_data);
			} else if ( sparse2->data != NULL ) {
				free(sparse2->data);
			}
			free(sparse2);
		}
		Mat_VarFree(pVariable);
	}

	void MatLabFile::Append( string_view pszFileName, string_view pszVariableName, const std::vector<string>& columnNames )
	{
		mat_t* pFile = Mat_Open( string(pszFileName).c_str(), MAT_ACC_RDWR );
		if ( !pFile )
			THROW( Exception("Could not create file '{}'"sv, pszFileName) );
		Append( pFile, pszVariableName, columnNames );
		Mat_Close( pFile );
	}
	void MatLabFile::Append( mat_t* pFile, string_view pszVariableName, const std::vector<string>& columnNames, size_t maxColumnNameLength/*=0*/ )
	{
		if( maxColumnNameLength==0 )
		{
			for( const auto& name : columnNames )
				maxColumnNameLength = max( maxColumnNameLength, name.length() );
		}
		size_t* columnDims = static_cast<size_t*>( calloc( 2,sizeof(size_t)) );
		columnDims[0] = columnNames.size();
		columnDims[1] = maxColumnNameLength; //columnNames.size();

		std::vector<char> columnNameArray( columnDims[0]*columnDims[1] );
		size_t iDataIndex = 0;
		for( size_t iColumn = 0;iColumn<columnDims[1]; ++iColumn )
		{
			for( size_t iRow=0; iRow<columnDims[0]; ++iRow )
			{
				const auto& columnName = columnNames[iRow];
				columnNameArray[iDataIndex++] = iColumn<columnName.length() ? columnName[iColumn] : ' ';
			}
		}
		MatLabVariable matvarColumns( pszVariableName, MAT_C_CHAR, matio_types::MAT_T_UINT8, 2, columnDims, (void*)columnNameArray.data(), 0 );
		matvarColumns.Write( pFile );
	}
/*	std::vector<string> MatLabFile::LoadColumns( string_view pszFileName, string_view pszVariableName )
	{
		//char* pszFileName = "C:\\Users\\duffyj\\Google Drive\\WorkShare\\Bosch\\Train\\numeric.mat";
		//char* pszVariableName = "x_columns";
		MatLabFile file( pszFileName );
		matvar_t* pVariable = Mat_VarRead( file._file, pszVariableName);
		const size_t rowCount = pVariable->dims[0];
		const size_t columnCount = pVariable->dims[1];
		std::vector<std::string> columnNames( rowCount );
		string_view pData = static_cast<string_view>( pVariable->data );
		for(int columnIndex=0, dataIndex=0; columnIndex<columnCount; ++columnIndex )
		{
			for(int rowIndex=0; rowIndex<rowCount; ++rowIndex )
			{
				string& value = columnNames[rowIndex];
				value += pData[dataIndex++];
			}
		}
		for( auto& columnName : columnNames )
		{
			StringUtilities::RTrim( columnName );
			clog << "columnName["<< columnName << "]" <<  endl;
		}
		return columnNames;
	}
	*/
}