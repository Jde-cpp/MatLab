#pragma once
#include <string>
#include <list>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Sparse>
//#include <matio.h>

#include "MatLabVariable.h"
//#include "../../../log/Logging.h"
//#include "../../../Exception.h"
//#include "../../../Stopwatch.h"

//#undef max
namespace std{
	class mutex;
	//class list;
	//class tuple;
}
/*
struct matvar_t;

struct _mat_t;
/ ** @brief Matlab MAT File information
* Contains information about a Matlab MAT file
* @ingroup MAT
/
typedef struct _mat_t mat_t;
*/

namespace Jde::IO::MatLab
{
//	using std::vector;
	class MatLabVariable;
#pragma warning( push )
#pragma warning( disable : 4251)
	class JDE_MATLAB_VISIBILITY MatLabFile
	{
	public:
		MatLabFile();
		//MatLabFile( const wchar_t* matlabFileName );
		MatLabFile( const fs::path& matlabFileName )noexcept(false);
		~MatLabFile();

		//matvar_t* operator[](std::size_t idx)       { return mVector[idx]; };
		/*matvar_t* operator[]( const std::string& variableName );*/
		MatLabVariable operator[]( string_view variableName );
		//MatLabVariable operator[]( string_view const variableName );

		void Combine( string_view path, string_view variableName, MatLabVariable v1, MatLabVariable v2 );

		template<typename T, typename CharType>
		static void SaveCsv( const std::basic_string<CharType>& csvFileName, const std::basic_string<CharType>& matFileName, const std::basic_string<CharType>& variabelName, const std::function<T(const std::basic_string<CharType>&)>& convertFunction, const std::function<T(const std::basic_string<CharType>&)>* pConvertFunction2=nullptr );

		template<typename T, typename CharType>
		static void LoadCsv( const std::basic_string<CharType>& csvFileName, const std::basic_string<CharType>& columnNamesToFetch, size_t maxLines );


		//template<typename T>
		//static void Save( const string_view& fileName, const string_view& variableName, const Eigen::SparseVector<T>& sparse );
		template<typename T>
		static void Save( string_view pszFileName, string_view variableName, const Eigen::SparseMatrix<T>& sparse, const std::vector<string>* pColumnNames=nullptr );
		template<typename T,int Rows,int Cols>
		static void Save( string_view pszFileName, string_view pszVariableName, const Eigen::Matrix<T, Rows, Cols>& matrix, const std::vector<string>* pColumnNames=nullptr );

		static void Append( string_view pszFileName, string_view pszVariableName, const std::vector<string>& columnNames );
		static void Append( mat_t* pFile, string_view pszVariableName, const std::vector<string>& columnNames, size_t maxColumnNameLength=0 );

		template<typename T, int Rows=-1, int Cols=-1>
		static void Append( string_view pszFileName, string_view pszVariableName, const Eigen::Matrix<T,Rows,Cols>& matrix );
		template<typename T>
		static void Append( string_view pszFileName, string_view pszVariableName, const Eigen::SparseMatrix<T>& sparse );

	private:
		static std::mutex* _lock;
		static void Logfunc( int log_level, char *message );
		static std::list< std::tuple<int,std::string> >* _logs;
		fs::path _fileName;

		mat_t* _file{NULL};
	};
#pragma region Save
	template<typename T,int Rows,int Cols>
	void MatLabFile::Save( string_view pszFileName, string_view pszVariableName, const Eigen::Matrix<T, Rows, Cols>& matrix, const std::vector<string>* pColumnNames/*=nullptr*/ )
	{
		mat_t* pFile = Mat_CreateVer( string(pszFileName).c_str(), nullptr, MAT_FT_MAT73 );
		if ( !pFile )
			THROW( Exception(fmt::format("Could not create file '{}'", pszFileName)) );

		MatLabVariable::Write( pFile, pszVariableName, matrix );
		if( pColumnNames )
			Append( pFile, string(pszVariableName)+"_columns", *pColumnNames );
		Mat_Close( pFile );
	}

	template<typename T>
	void MatLabFile::Save( string_view fileName, string_view variableName, const Eigen::SparseMatrix<T>& sparse, const std::vector<string>* pColumnNames/*=nullptr*/ )
	{
		//Save( boost::locale::conv::utf_to_utf<wchar_t,char>(fileName), boost::locale::conv::utf_to_utf<wchar_t,char>(variableName), sparse );
		mat_t* pFile = Mat_CreateVer( string(fileName).c_str(), NULL, MAT_FT_MAT73 );
		if ( !pFile )
			THROW( Exception(fmt::format("Could not create file '{}'", fileName)) );
		DBG( "Saving:  {} [{}]", fileName, variableName );
		MatLabVariable::Write( pFile, variableName, sparse );
		if( pColumnNames )
			Append( pFile, string(variableName)+"_columns", *pColumnNames );

		Mat_Close( pFile );

	}

	//template<typename T>
	//void MatLabFile::Save( const string_view& fileName, const string_view& variableName, const Eigen::SparseMatrix<T>& sparse )
	//{
	//	mat_t* pFile = Mat_CreateVer( boost::locale::conv::utf_to_utf<char,wchar_t>(fileName).c_str(), NULL, MAT_FT_MAT73 );
	//	if ( !pFile )
	//		THROW( Exception(boost::wformat{L"Could not create file '%1%'"} % fileName) );
	//	BOOST_LOG_TRIVIAL(debug) << "Saving:  " << fileName << "[" << variableName << "]" << std::endl;
	//	MatLabVariable::Write( pFile, variableName, sparse );

	//	Mat_Close( pFile );
	//}
	template<typename T, int Rows, int Cols>
	void MatLabFile::Append( string_view pszFileName, string_view pszVariableName, const Eigen::Matrix<T,Rows,Cols>& matrix )
	{
		mat_t* pFile = Mat_Open( string(pszFileName).c_str(), MAT_ACC_RDWR );
		if ( !pFile )
			THROW( Exception(fmt::format("Could not create file '{}'", pszFileName)) );
		TRACE( "Saving:  {} [{}]"sv, pszFileName, pszVariableName );
		MatLabVariable::Write<T,Rows,Cols>( pFile, string(pszVariableName).c_str(), matrix );

		Mat_Close( pFile );
	}

	template<typename T>
	void MatLabFile::Append( string_view  pszFileName, string_view  pszVariableName, const Eigen::SparseMatrix<T>& sparse )
	{
		mat_t* pFile = Mat_Open( string(pszFileName).c_str(), MAT_ACC_RDWR );
		if ( !pFile )
			THROW( Exception(fmt::format("Could not create file '{}'",pszFileName)) );
		TRACE( "Saving:  {} [{}]", pszFileName, pszVariableName );
		MatLabVariable::Write( pFile, pszVariableName, sparse );

		Mat_Close( pFile );
	}

#pragma endregion
#pragma region SaveCsv

	template<typename T, typename CharType>
	void MatLabFile::SaveCsv( const std::basic_string<CharType>& csvFileName, const std::basic_string<CharType>& matFileName, const std::basic_string<CharType>& variableName, const std::function<T(const std::basic_string<CharType>&)>& convertFunction, const std::function<T(const std::basic_string<CharType>&)>* pConvertFunction2 )
	{
		Stopwatch sw( StopwatchTypes::ReadFile, variableName.c_str(), true );
		std::vector<string> columnNames;
		size_t maxColumnNameLength=0;
		constexpr size_t maxLines = 0;
		vector<int> columnCounts;
		size_t lineCount=0;
		auto getStats = [&sw, &lineCount,&columnCounts,&columnNames,&maxColumnNameLength/*,&columnNameLength*/](const basic_string<CharType>& line)mutable
		{
			const auto columnCount =  count(line.begin(), line.end(), ',');
			list<std::basic_string<CharType>> tokens = StringUtilities::Split<CharType>(line);
			if( lineCount==0 )
			{
				columnCounts.resize( columnCount+1 );
				for( const auto& token : tokens )
				{
					maxColumnNameLength = std::max<size_t>( maxColumnNameLength, token.length() );
					columnNames.push_back( token );
					//columnNameLength+=token.length()+1;
				}
			}
			else
			{
				if( tokens.size()!=columnCounts.size() )
					THROW( Exception(fmt::format("Column counts don't add up for line '{}' actual:  '{}' expected:  '{}'", lineCount, tokens.size(), columnCounts.size())) );

				auto columnIndex = 0;
				for( auto pToken = tokens.begin(); pToken!=tokens.end(); ++columnIndex, ++pToken )
				{
					if( !pToken->empty() )
						++columnCounts[columnIndex];
				}
			}
			if( ++lineCount%50000==0 )
				sw.Progress( lineCount );
		};
		if( maxLines==0 )
			IO::File::ForEachLine<CharType>( csvFileName, getStats );
		else
			IO::File::ForEachLine<CharType>( csvFileName, getStats, maxLines );

		//sw.Finish();
		Stopwatch sw2( StopwatchTypes::ReadFile, variableName.c_str(), true );
		int valueCount = 0;
		for( auto columnIndex=0; columnIndex<columnCounts.size(); ++columnIndex )
			valueCount+=columnCounts[columnIndex];

		T* values = static_cast<T*>( calloc(valueCount, sizeof(T)) );
		int lineIndex=0;
		int* ir = static_cast<int*>( calloc( valueCount, sizeof(int)) );

		//size_t maxDataIndex = numeric_limits<size_t>::min();
		//size_t minDataIndex = numeric_limits<size_t>::max();
		vector<int> columnValueCount(columnCounts.size());
		auto getValues = /*[&]*/[&ir,&values, &lineIndex,&columnCounts,&columnValueCount, &convertFunction, &sw2, &lineCount, pConvertFunction2](const std::basic_string<CharType>& line)mutable
		{
			if( lineIndex!=0 )
			{
				list<std::basic_string<CharType>> tokens = StringUtilities::Split(line);
				auto columnIndex = 0;
				for( auto pToken = tokens.begin(); pToken!=tokens.end(); ++columnIndex, ++pToken )
				{
					if( pToken->empty() )
						continue;
					/*double value = NumberUtilities::TryToDouble(*pToken);
					if( std::isnan(value) )
						continue;
						*/
					T value = columnIndex==0 && pConvertFunction2 ? (*pConvertFunction2)(*pToken) : convertFunction(*pToken);
					size_t dataIndex = 0;
					for( size_t previousIndex=0; previousIndex<columnIndex;++previousIndex )
						dataIndex+=columnCounts[previousIndex];
					dataIndex+=columnValueCount[columnIndex];
					ir[dataIndex] = lineIndex-1;
					values[dataIndex] = value;
					//maxDataIndex = max(maxDataIndex, dataIndex);
					//minDataIndex = min(minDataIndex, dataIndex);

					++columnValueCount[columnIndex];
				}
			}
			if( ++lineIndex%50000==0 )
				sw2.Progress( lineIndex, lineCount );
		};
		if( maxLines==0 )
			IO::File::ForEachLine<CharType>( csvFileName, getValues );
		else
			IO::File::ForEachLine<CharType>( csvFileName, getValues, maxLines );

		sw2.Finish();
		Stopwatch sw3( StopwatchTypes::WriteFile, variableName.c_str(), true );
		int* jc = static_cast<int*>( calloc((columnCounts.size()+1),sizeof(int)) );
		jc[0] = 0;
		for( int columnIndex=0; columnIndex<columnCounts.size(); ++columnIndex )
			jc[columnIndex+1] = jc[columnIndex]+columnCounts[columnIndex];

		mat_sparse_t* sparse_data = static_cast<mat_sparse_t*>( calloc(1, sizeof(mat_sparse_t)) );// = new(mat_sparse_t*)calloc( 1, sizeof(*sparse_data) );
		sparse_data->nir = sparse_data->ndata = sparse_data->nzmax = valueCount;
		sparse_data->ir = ir;
		sparse_data->jc = jc;
		sparse_data->njc = int( columnCounts.size()+1 );
		sparse_data->data = values;

		DBG( "dataCount:  {}", sparse_data->ndata );
		DBG( "ColumnCount:  {}", columnCounts.size() );
//		clog << "zeroCount:  " << zeroCount << endl;
//		clog << "maxDataIndex:  " << maxDataIndex << endl;
//		clog << "minDataIndex:  " << minDataIndex << endl;
//		clog << "ColumnCount:  " << columnCounts.size() << endl;
		//clog << "RowCount:  " << lineCount-1 << endl;

		size_t* dims = static_cast<size_t*>( calloc( 2,sizeof(size_t)) );
		dims[0] = lineCount-1;
		dims[1] = columnCounts.size();

		const matio_types dataType = MatLabVariable::GetMatioType( typeid(T) );
		MatLabVariable matvar( variableName, MAT_C_SPARSE, dataType, 2, dims, static_cast<void*>(sparse_data), 0 );

		mat_t* pFile = Mat_CreateVer( matFileName.c_str(), NULL, MAT_FT_MAT73/*MAT_FT_DEFAULT*/ );
		if ( !pFile )
			THROW( Exception(fmt::format("Could not create file '{1}'", matFileName)) );

		//Mat_VarWrite( pFile, _pVariable, MAT_COMPRESSION_ZLIB );
		matvar.Write( pFile );

		Append( pFile, (variableName+"_columns").c_str(), columnNames, maxColumnNameLength );
		Mat_Close( pFile );
	}
#pragma endregion
#pragma warning( pop )
}