#pragma once
#pragma warning(disable :4127)
#include <Eigen/Dense>
#include <Eigen/Sparse>
#pragma warning(default :4127)
#include <matio.h>
#include "Exports.h"
#include "../../Framework/source/TypeDefs.h"
#include "../../Framework/source/Exception.h"
#include "../../Framework/source/Stopwatch.h"

//struct matvar_t;
namespace Jde::IO::MatLab
{
	class MatLabVariable
	{
	public:
		friend class MatLabFile;

		JDE_MATLAB_VISIBILITY MatLabVariable( string_view name, const Eigen::MatrixXd& matrix );
		JDE_MATLAB_VISIBILITY MatLabVariable( string_view name, const Eigen::VectorXd& vector );
		JDE_MATLAB_VISIBILITY MatLabVariable( string_view name, const std::shared_ptr<Eigen::SparseMatrix<double>> matrix );
		//JDE_MATLAB_VISIBILITY MatLabVariable( const std::wstring& name, matio_classes classType, matio_types data_type, int rank, size_t *dims, void *data, int opt );
		JDE_MATLAB_VISIBILITY MatLabVariable( string_view name, matio_classes classType, matio_types data_type, int rank, size_t *dims, void *data, int opt );
		JDE_MATLAB_VISIBILITY MatLabVariable( matvar_t* variable );
		JDE_MATLAB_VISIBILITY MatLabVariable( const MatLabVariable& copyFrom );
		JDE_MATLAB_VISIBILITY ~MatLabVariable();

		//JDE_MATLAB_VISIBILITY sp<Eigen::MatrixXd> ToMatrix()const noexcept(false);
		//JDE_MATLAB_VISIBILITY sp<Eigen::VectorXd> ToVector()const noexcept(false);

		JDE_MATLAB_VISIBILITY std::vector<std::string> ToColumnNames();
		template<typename T,int Rows=-1,int Cols=-1>
		std::unique_ptr<Eigen::Matrix<T,Rows,Cols>> ToEigen()const;
		template<int Rows,int Cols>
		std::vector<std::unique_ptr<Eigen::Matrix<double,Rows,Cols>>> ToEigenVector()const;

		template<typename Scaler>
		static void ToSparseReserveMatrix( Eigen::SparseMatrix<Scaler>* result, const Eigen::VectorXi& columnCounts, int rowCount );
		template<typename Scaler>
		static void ToSparseReserveVector( Eigen::SparseVector<Scaler>* result, const Eigen::VectorXi& columnCounts, int rowCount );
		template<typename Result,typename Scaler,typename FileType>
		std::unique_ptr<Result> ToSparse( const Scaler* zeroTest, const Scaler* zeroValue, const std::function<void(Result*,const Eigen::VectorXi&,int)>& reserveFunction, size_t rowCount=std::numeric_limits<size_t>::max() )const;

		JDE_MATLAB_VISIBILITY double operator()(size_t rowIndex, size_t columnIndex) const;

		JDE_MATLAB_VISIBILITY void Write( mat_t* pFile )const;

		template<typename T>
		static void Write( mat_t* pFile, string_view variableName, const Eigen::SparseMatrix<T>& sparse );

		template<typename T>
		static void Write( mat_t* pFile, string_view variableName, const Eigen::SparseVector<T>& sparse );

		JDE_MATLAB_VISIBILITY static matio_types GetMatioType( const std::type_info& typeInfo );
		JDE_MATLAB_VISIBILITY static matio_classes GetMatioClass( const std::type_info& typeInfo );
		JDE_MATLAB_VISIBILITY const std::string& GetName()const noexcept{ return _name;}
		JDE_MATLAB_VISIBILITY const size_t rows()const noexcept(false){return GetVariable()->dims[0];}
		JDE_MATLAB_VISIBILITY const size_t cols()const noexcept(false){return GetVariable()->dims[1];}
		matvar_t* GetVariable(){ assert(_pVariable!=nullptr); return _pVariable;}
		const matvar_t* GetVariable()const{ assert(_pVariable!=nullptr); return _pVariable;}

		JDE_MATLAB_VISIBILITY const std::vector<std::pair<size_t,size_t>> ColumnCounts()const noexcept(false);

		JDE_MATLAB_VISIBILITY const mat_sparse_t* GetSparsePointer()const noexcept(false);


		template<typename ScalerFrom, typename ScalerTo>
		void ForEachValue( std::function<void( size_t rowIndex,size_t columnIndex, ScalerTo value )> func, Stopwatch* sw  )const noexcept(false);
	private:
		template<typename T,int Rows,int Cols>
		static void Write( mat_t* pFile, string_view  pszVariableName, const Eigen::Matrix<T, Rows, Cols>& matrix );

		mutable matvar_t* _pVariable{NULL};
		std::string _name;
		enum MatrixType{Data,Matrix,Vector, Sparse};
		MatrixType _matrixType;
		union
		{
			const Eigen::MatrixXd* _matrix;
			const Eigen::VectorXd* _vector;
			std::shared_ptr<Eigen::SparseMatrix<double>> _sparse;
		};
	};
#pragma region Write
	template<typename T,int Rows,int Cols>
	void MatLabVariable::Write( mat_t* pFile, string_view  pszVariableName, const Eigen::Matrix<T, Rows, Cols>& matrix )
	{
		const size_t rowCount = matrix.rows();
		const size_t columnCount = matrix.cols();
		size_t* dims = static_cast<size_t*>( malloc( 2*sizeof(size_t)) );
		dims[0] = rowCount;
		dims[1] = columnCount;
		//DBG( "rowCount={}, columnCount={}", rowCount, columnCount );
		T* values = static_cast<T*>( calloc(rowCount*columnCount, sizeof(T)) );
		//matlab file is column major & eigen defaults to column major
		for( uint columnIndex=0; columnIndex<columnCount; ++columnIndex )
		{
			for( uint rowIndex=0; rowIndex<rowCount; ++rowIndex )
				values[columnIndex*rowCount+rowIndex] = matrix(rowIndex, columnIndex);
		}
		matvar_t* pVariable = Mat_VarCreate( string(pszVariableName).c_str(), GetMatioClass(typeid(T)),GetMatioType(typeid(T)),2,dims,values,0 );
		Mat_VarWrite( pFile, pVariable, MAT_COMPRESSION_ZLIB/*MAT_COMPRESSION_NONE*/ );
		Mat_VarFree(pVariable);
		std::free( values );
	}
	template<typename T>
	void MatLabVariable::Write( mat_t* pFile, string_view variableName, const Eigen::SparseVector<T>& sparse )
	{
		Stopwatch sw( variableName );
		const auto columnCount = sparse.cols();
		auto jc = static_cast<mat_uint32_t*>( calloc(2, sizeof(mat_uint32_t)) );
		const auto valueCount = sparse.nonZeros();
		jc[1] = valueCount;
		//std::clog << "valueCount:  " << valueCount << ".  rows:  " << sparse.rows() << endl;
		T* values = static_cast<T*>( calloc(valueCount,sizeof(T)) );
		auto ir = static_cast<mat_uint32_t*>( calloc( valueCount,sizeof(mat_uint32_t)) );
		if( ir==nullptr )
			THROW( Exception( "out of memory" ) );
		auto dataIndex = 0;
		/*TODO
		for( Eigen::SparseVector<T>::InnerIterator it(sparse,0); it; ++it )
		{
			const auto rowIndex = it.row();
			ir[dataIndex] = rowIndex;
			const T value = it.value();
			values[dataIndex++] = value==std::numeric_limits<T>::epsilon() ? 0 : value;
		}*/
		//std::clog << "dataIndex:  " << dataIndex << endl;
		assert( dataIndex==valueCount );
		mat_sparse_t* sparse_data = static_cast<mat_sparse_t*>( calloc(1, sizeof(mat_sparse_t)) );
		sparse_data->nir = sparse_data->ndata = sparse_data->nzmax = valueCount;
		sparse_data->ir = ir;

		sparse_data->jc = jc;
		sparse_data->njc = columnCount+1;
		sparse_data->data = values;

		size_t* dims = static_cast<size_t*>( malloc( 2*sizeof(size_t)) );
		dims[0] = sparse.rows();
		dims[1] = columnCount;
		auto pVariable = Mat_VarCreate( string(variableName).c_str(), MAT_C_SPARSE,GetMatioType( typeid(T) ),2,dims,static_cast<void*>(sparse_data), MAT_F_DONT_COPY_DATA );
		//try
		//{
		Mat_VarWrite( pFile, pVariable, MAT_COMPRESSION_ZLIB );
		//Logging::ProcessInfo::LogMemoryUsage( wstring("before Mat_VarFree") );
		//pVariable->mem_conserve = 0;
		{
			mat_sparse_t* sparse2 = static_cast<mat_sparse_t*>(pVariable->data );
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

	template<typename T>
	void MatLabVariable::Write( mat_t* pFile, string_view variableName, const Eigen::SparseMatrix<T>& sparse )
	{
		Stopwatch sw( fmt::format("StopwatchTypes::WriteFile - {}", variableName), true );
		const int64_t columnCount = sparse.cols();
		mat_uint32_t* jc = static_cast<mat_uint32_t*>( calloc(columnCount+1, sizeof(mat_uint32_t)) );
		//for( int columnIndex=0; columnIndex<columnCount; ++columnIndex )
		//			jc[columnIndex+1] = jc[columnIndex]+sparse.innerNonZeroPtr()[columnIndex];
		const int64_t valueCount = sparse.nonZeros();
		DBG( "valueCount:  {}"sv, valueCount );
		T* values = static_cast<T*>( calloc(valueCount,sizeof(T)) );
		mat_uint32_t* ir = static_cast<mat_uint32_t*>( calloc( valueCount,sizeof(mat_uint32_t)) );
		if( ir==nullptr )
			THROW( Exception( "out of memory" ) );
		auto dataIndex = 0;
		for( auto columnIndex=0; columnIndex<sparse.outerSize(); ++columnIndex )
		{
			jc[columnIndex+1] = jc[columnIndex];
			/* TODOI
			for( Eigen::SparseMatrix<T>::InnerIterator it(sparse,columnIndex); it; ++it )
			{
				++jc[columnIndex+1];
				const auto rowIndex = it.row();
				ir[dataIndex] = rowIndex;
				const T value = it.value();
				values[dataIndex++] = value==std::numeric_limits<T>::epsilon() ? 0 : value;
			}
			*/
		}
		assert( dataIndex==valueCount );
		DBG( "dataIndex:  {}."sv, dataIndex );
		mat_sparse_t* sparse_data = static_cast<mat_sparse_t*>( calloc(1, sizeof(mat_sparse_t)) );
		sparse_data->nir = sparse_data->ndata = sparse_data->nzmax = static_cast<int>( valueCount );
		sparse_data->ir = ir;

		sparse_data->jc = jc;
		sparse_data->njc = static_cast<int>( columnCount+1 );
		sparse_data->data = values;

		size_t* dims = static_cast<size_t*>( malloc( 2*sizeof(size_t)) );
		dims[0] = sparse.rows();
		dims[1] = columnCount;
		auto pVariable = Mat_VarCreate( string(variableName).c_str(), MAT_C_SPARSE,GetMatioType( typeid(T) ),2,dims,static_cast<void*>(sparse_data), MAT_F_DONT_COPY_DATA );
		//try
		//{
		Mat_VarWrite( pFile, pVariable, MAT_COMPRESSION_ZLIB );
		//Logging::ProcessInfo::LogMemoryUsage( wstring("before Mat_VarFree") );
		//pVariable->mem_conserve = 0;
		{
			mat_sparse_t *sparse2 = static_cast<mat_sparse_t*>(pVariable->data);
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
#pragma endregion
#pragma region ToSparse
	template<typename Scaler>
	void MatLabVariable::ToSparseReserveMatrix( Eigen::SparseMatrix<Scaler>* result, const Eigen::VectorXi& columnCounts, int rowCount )
	{
		result->reserve( columnCounts );
	}

	template<typename Scaler>
	void MatLabVariable::ToSparseReserveVector( Eigen::SparseVector<Scaler>* result, const Eigen::VectorXi& columnCounts, int rowCount )
	{
		const int size = rowCount;
		result->resize( size );//needs to change m_size
		result->reserve( size );
	}

	template<typename Result,typename Scaler,typename FileType>
	std::unique_ptr<Result> MatLabVariable::ToSparse(const Scaler* zeroTest, const Scaler* zeroValue, const std::function<void(Result*,const Eigen::VectorXi&,int)>& reserveFunction, size_t maxRowCount )const
	{
		Stopwatch sw( _name );
		auto pVariable = GetVariable();
		//const size_t stride = Mat_SizeOf( pVariable->data_type );
		if( pVariable->class_type!=MAT_C_SPARSE )
			THROW( Exception("class_type!=MAT_C_SPARSE") );

		const mat_sparse_t* sparse =  GetSparsePointer();
		std::ostringstream os;
		const auto columnCount = sparse->njc-1;
		const auto fileRowCount= int(pVariable->dims[0]);
		const auto rowCount= std::min<size_t>( fileRowCount, maxRowCount );
		Eigen::VectorXi columnCounts( columnCount );
		if( rowCount<fileRowCount )
		{
			columnCounts = Eigen::VectorXi::Zero( columnCount );
			for( auto jcIndex=1; jcIndex<sparse->njc; ++jcIndex )
			{
				const auto columnIndex = jcIndex-1;
				for( auto dataIndex = sparse->jc[columnIndex]; dataIndex<sparse->jc[jcIndex]; ++dataIndex )
				{
					const auto rowIndex = sparse->ir[dataIndex];
					if( rowIndex>=rowCount )
						break;
					++columnCounts( columnIndex );
				}
			}
		}
		else
		{
			for( auto jcIndex=1, startIndex=0; jcIndex<sparse->njc; ++jcIndex )
			{
				const auto currentIndex = sparse->jc[jcIndex];
				columnCounts( jcIndex-1 ) = currentIndex-startIndex;
				startIndex = currentIndex;
			}
		}

		Result* matrix = new Result( int(rowCount), columnCount );
		reserveFunction( matrix, columnCounts, int(rowCount) );
		FileType* data = static_cast<FileType*>(sparse->data);
		const bool testZero = zeroTest!=nullptr;
		const size_t progressIndex = std::max(10000u,sparse->nir/10);
		for( auto jcIndex=1; jcIndex<sparse->njc; ++jcIndex )
		{
			const auto columnIndex = jcIndex-1;
			for( auto dataIndex = sparse->jc[jcIndex-1]; dataIndex<sparse->jc[jcIndex]; ++dataIndex )
			{
				const auto rowIndex = sparse->ir[dataIndex];
				if( rowIndex>=rowCount )
					break;
				Scaler value =  Scaler( data[dataIndex] );

				if( testZero && value==*zeroTest )
					value = *zeroValue;//std::numeric_limits<T>::epsilon()
				matrix->insert(rowIndex,columnIndex) = value;
				if( (dataIndex+1)%progressIndex==0 )
					sw.Progress( dataIndex, sparse->nir );
			}
		}

		return std::unique_ptr<Result>( matrix );
	}
#pragma endregion
#pragma region ForEachValue
	template<typename ScalerFrom, typename ScalerTo>
	void MatLabVariable::ForEachValue( std::function<void( size_t rowIndex,size_t columnIndex, ScalerTo value )> func, Stopwatch* sw  )const
	{
		if( _pVariable->class_type==MAT_C_SPARSE )
		{
			const auto sparse = GetSparsePointer();
			const size_t progressIndex = std::max( 10000u, sparse->nir/10 );
			const ScalerFrom* pData = static_cast<ScalerFrom*>( sparse->data );
			for( auto jcIndex=1; jcIndex<sparse->njc; ++jcIndex )
			{
				const auto columnIndex = jcIndex-1;
				for( auto dataIndex = sparse->jc[jcIndex-1]; dataIndex<sparse->jc[jcIndex]; ++dataIndex )
				{
					const auto rowIndex = sparse->ir[dataIndex];
					ScalerTo value =  ScalerTo( pData[dataIndex] );
					func( rowIndex, columnIndex, value );
					if( sw && (dataIndex+1)%progressIndex==0 )
						sw->Progress( dataIndex, sparse->nir );
				}
			}
		}
		else
		{
			auto pVariable = GetVariable();
			//const size_t dims =  pVariable->dims[2];
			const size_t rowCount =  rows();
			const size_t columnCount =  cols();
			const size_t stride = Mat_SizeOf( _pVariable->data_type );
			const char* pData = static_cast<const char*>( pVariable->data );
			for( size_t rowIndex = 0; rowIndex < rowCount; ++rowIndex )
			{
				for ( size_t columnIndex = 0; columnIndex < columnCount; ++columnIndex )
				{
					const size_t index = rowCount*columnIndex+rowIndex;
					ScalerFrom* value = (ScalerFrom*)( pData+index*stride );
					func( rowIndex, columnIndex, ScalerTo(*value) );
				}
			}
		}
	}
#pragma endregion
#pragma region ToEigen
	template<typename T,int Rows,int Cols>
	std::unique_ptr<Eigen::Matrix<T,Rows,Cols>> MatLabVariable::ToEigen()const
	{
		//MatLab::MatLabVariable variable = file[variableName];
		if( _pVariable==NULL )
			throw Jde::Exception( "Variable not initialized" );

		size_t stride = Mat_SizeOf( _pVariable->data_type );

	//	if(  _pVariable->class_type!=MAT_C_DOUBLE )
	//		throw Exception( "class_type!=MAT_C_DOUBLE" );

		const uint rowCount =  _pVariable->dims[0];
		const uint columnCount =  _pVariable->dims[1];
		auto pMatrix = new Eigen::Matrix<T,Rows,Cols>( rowCount, columnCount );
		for( uint i = 0; i < rowCount; ++i )
		{
			for ( uint j = 0; j < columnCount; ++j )
			{
				const auto index = rowCount*j+i;
				const char* data = static_cast<const char*>(_pVariable->data);
				const T* value = (T*)( data+index*stride );
				(*pMatrix)(i,j) = *value;
			}
		}
		return std::unique_ptr<Eigen::Matrix<T,Rows,Cols>>( pMatrix );
	}
#pragma endregion
#pragma region ToEigenVector
	template<int Rows,int Cols>
	std::vector<std::unique_ptr<Eigen::Matrix<double,Rows,Cols>>> MatLabVariable::ToEigenVector()const
	{
		if( _pVariable==NULL )
			throw Jde::Exception( "Variable not initialized" );

		size_t stride = Mat_SizeOf( _pVariable->data_type );

		if( _pVariable->class_type!=MAT_C_DOUBLE )
			throw Exception( "class_type!=MAT_C_DOUBLE" );

		const size_t rowCount =  _pVariable->dims[0];
		const size_t columnCount =  _pVariable->dims[1];
		const size_t dims =  _pVariable->dims[2];
		const char* data = static_cast<const char*>(_pVariable->data);

		std::vector<std::unique_ptr<Eigen::Matrix<double,Rows,Cols>>> matrixes = std::vector<std::unique_ptr<Eigen::Matrix<double,Rows,Cols>>>(dims);
		for ( size_t dimIndex = 0; dimIndex < dims; ++dimIndex )
			matrixes[dimIndex] = std::unique_ptr<Eigen::Matrix<double,Rows,Cols>>( new Eigen::Matrix<double,Rows,Cols>(rowCount,columnCount) );

		for( size_t rowIndex = 0; rowIndex < rowCount; ++rowIndex )
		{
			for ( size_t columnIndex = 0; columnIndex < columnCount; ++columnIndex )
			{
				for ( size_t dimIndex = 0; dimIndex < dims; ++dimIndex )
				{
					size_t index = dimIndex*(rowCount*columnCount)+rowCount*columnIndex+rowIndex;
					double* value = (double*)( data+index*stride );
					(*matrixes[dimIndex])(rowIndex,columnIndex) = *value;
				}
			}
		}
		return matrixes;
	}
#pragma endregion
}