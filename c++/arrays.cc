#include <iostream>
#include <new>

using namespace std;

//------------------------------------------------------------------------------
void ErrorMessage(const char* message)
{
   cout << message << endl;
   exit(1);
}

//------------------------------------------------------------------------------
// Allocate a 1D array
//------------------------------------------------------------------------------
template <class Type> 
void Allocate(Type *&m, int d1)
{
	m = new (nothrow) Type [d1];
	if (m == 0)
		ErrorMessage("Error: Memory can not be allocated.");

	for (int i = 0; i < d1; i++)
		m[i] = (Type)(0.0);
}

//------------------------------------------------------------------------------
// Deallocate a 1D array
//------------------------------------------------------------------------------
template <class Type> 
void Deallocate(Type *m, int d1)
{
	delete [] m;
	m = NULL;
}

//------------------------------------------------------------------------------
// Allocate a 2D array
//------------------------------------------------------------------------------
template <class Type> 
void Allocate(Type **&m, int d1, int d2)
{
	m = new (nothrow) Type* [d1];
	
	if (m == 0)
		ErrorMessage("Error: Memory can not be allocated.");
	
	for (int i = 0; i < d1; i++)
	{
      m[i] = new (nothrow) Type [d2];
		
		if (m[i] == 0)
			ErrorMessage("Error: Memory can not be allocated.");

		for (int j = 0; j < d2; j++)
		{
			m[i][j] = (Type)(0.0);
      }
	}
}

//------------------------------------------------------------------------------
// Deallocate a 2D array
//------------------------------------------------------------------------------
template <class Type> 
void Deallocate(Type **m, int d1, int d2)
{
	for (int i = 0; i < d1; i++)
   {
			delete [] m[i];
	}
	
	delete [] m;
	m = NULL;
}

//------------------------------------------------------------------------------
// Allocate a 3D array
//------------------------------------------------------------------------------
template <class Type> 
void Allocate(Type ***&m, int d1, int d2, int d3)
{
	m = new (nothrow) Type** [d1];
	
	if (m == 0)
		ErrorMessage("Error: Memory can not be allocated.");
	
	for (int i = 0; i < d1; i++)
	{
      m[i] = new (nothrow) Type* [d2];
		
		if (m[i] == 0)
			ErrorMessage("Error: Memory can not be allocated.");

		for (int j = 0; j < d2; j++)
		{
			m[i][j] = new (nothrow) Type [d3];
			
			if (m[i][j] == 0)
			   ErrorMessage("Error: Memory can not be allocated.");
			
			for (int k = 0; k < d3; k++)
				m[i][j][k] = (Type)(0.0);
		}
	}
}

//------------------------------------------------------------------------------
// Deallocate a 3D array
//------------------------------------------------------------------------------
template <class Type> 
void Deallocate(Type ***m, int d1, int d2, int d3)
{
	for (int i = 0; i < d1; i++)
   {
		for (int j = 0; j < d2; j++)
			delete [] m[i][j];
		delete [] m[i];
	}
	
	delete [] m;
	m = NULL;
}

//------------------------------------------------------------------------------
int main()
{
   int n1 = 10;
   int n2 = 20;
   int n3 = 30;

   // 1d array
   double* a;
   Allocate(a, n1);

   // 2d array
   double** b;
   Allocate(b, n1, n2);

   // 3d array
   double*** c;
   Allocate(c, n1, n2, n3);

   // You must deallocate, otherwise memory will build up.
   // And you need to remember the size of each array.
   Deallocate(a, n1);
   Deallocate(b, n1, n2);
   Deallocate(c, n1, n2, n3);

	return 0;
}
