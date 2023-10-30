#include "plot.h"
//-------Global variables relevant for graphics-----------------------

  extern float gfa_xmold[];	//MAXDIM

  char  gca_sbuft[80];        // Stringbuffer

  float gf_xgraph,gf_ygraph,gf_ygo,gf_xgo;
  float gf_x_low, gf_x_high, gf_y_low, gf_y_high;
  float gf_x_l, gf_x_h, gf_y_l, gf_y_h;

  static int gi_zoom_flag = 0;
  static float gf_x_low_temp, gf_x_high_temp, gf_y_low_temp, gf_y_high_temp;
  static float gf_x_move, gf_y_move;
  static float gf_x_low_sav, gf_x_high_sav, gf_y_low_sav, gf_y_high_sav;


//-----Initialization of graphics-------------------------------------

void graphics_init(void)
/**C*F****************************************************************
**                                                                  
** Function       :void graphics_init()                                        
**                                                                  
** Author         :Rainer Storn                                     
**                                                                  
** Descripton     :Defines the coordinate system.                 
**                                                                  
** Functions      :axis()                                                
**                                                                  
** Globals        :(see code)                                                
**                                                                  
** Parameters     :-   
**                                                                  
** Preconditions  :-                     
**                                                                  
** Postconditions :-       
**
** Return Value   :-                                                                                  
**                                                                  
***C*F*E*************************************************************/
{
  // define coordinate system
  gf_x_low  = -1.5;
  gf_x_high =  1.5;
  gf_y_low  = -5;
  gf_y_high = 10;

  // save the coordinate system to allow zooming
  gf_x_low_sav  = gf_x_low;
  gf_x_high_sav = gf_x_high;
  gf_y_low_sav  = gf_y_low;
  gf_y_high_sav = gf_y_high;

  // define axes
  axis(gf_x_low,gf_x_high,gf_y_low,gf_y_high); 

  gca_sbuft[0] = '\0';            // Initialize string buffer
}

void draw_graph(float fa_params[], int i_D, char color)
/**C*F****************************************************************
**                                                                  
** Function       :void draw_graph()                                        
**                                                                  
** Author         :Rainer Storn                                     
**                                                                  
** Description    :-                 
**                                                                  
** Functions      :fline()                                                
**                                                                  
** Globals        :(see code)                                                
**                                                                  
** Parameters     :fa_params[]   (I)    parameter vector
**                 i_D           (I)    size of parameter vector
**                 color         (I)    character variable defining the color   
**                                                                  
** Preconditions  :graphics_init() must have been called.                     
**                                                                  
** Postconditions :- 
** 
** Return Value   :-                                                                                  
**                                                                                                                                    
***C*F*E*************************************************************/
{
   int i,j;

//------first point--------------------------------
		gf_xgo = -1.2;
		gf_ygraph = 0.0;
		for (j=0; j<i_D-1; j++)
		{
		  gf_ygraph = (gf_ygraph + fa_params[j])*gf_xgo;
		}
		gf_ygraph+=fa_params[i_D-1];

        if (gf_ygraph < gf_y_low)  gf_ygraph = gf_y_low;
        if (gf_ygraph > gf_y_high) gf_ygraph = gf_y_high;
		gf_ygo = gf_ygraph;

//------Compute remaining points.------------------

		for (i=1; i<=120; i++)
		{
		  gf_ygraph = 0.0;
		  gf_xgraph = -1.2 + (float)i/50;
		  for (j=0; j<i_D-1; j++)
		  {
		  gf_ygraph = (gf_ygraph + fa_params[j])*gf_xgraph;
		  }
		  gf_ygraph+=fa_params[i_D-1];

		  if (gf_xgraph < gf_x_low)  gf_xgraph = gf_x_low;
		  if (gf_xgraph > gf_x_high) gf_xgraph = gf_x_high;
		  if (gf_ygraph < gf_y_low)  gf_ygraph = gf_y_low;
		  if (gf_ygraph > gf_y_high) gf_ygraph = gf_y_high;

		  fline(gf_xgo,gf_ygo,gf_xgraph,gf_ygraph,color);
		  gf_ygo = gf_ygraph;
		  gf_xgo = gf_xgraph;
		}
}

void update_graphics(float best[], int i_D, float fa_bound[], long l_nfeval, int i_gen, float f_emin, int i_strategy, int gi_genmax)
/**C*F****************************************************************
**                                                                  
** Function       :void update_graphics()                                        
**                                                                  
** Author         :Rainer Storn                                     
**                                                                  
** Description    :Custom program which updates the graphics part of the
**                 differential evolution optimization.                 
**                                                                  
** Functions      :xflt(), yflt(), frect(), Cls(), axis(), draw_graph(), box(),
**                 grid(), fline(), sprintf(), myprint().                                                
**                                                                  
** Globals        :(see code)  
**                                                                  
** Parameters     :best[]        (I)    parameter vector
**                 i_D           (I)    dimension of the parameter vector  
**                 fa_bound[]    (I)    array defining a tolerance scheme for the current example
**                 l_nfeval      (I)    current number of acumulated function evaluations                                          
**                 i_gen         (I)    current number of accumulated generations
**                 f_emin        (I)    current best objective function value
**                 i_strategy    (I)    DE-strategy used (coded as a number)   
**                                                                  
** Preconditions  :graphics_init() must have been called.                     
**                                                                  
** Postconditions :-                                             
**                  
** Return Value   :-                                                                                  
**                                                                                                                   
***C*F*E*************************************************************/
{
	int j;

//======Zoom control================================================

	if ((MouseLp == TRUE) && (MouseL == TRUE) && (gi_zoom_flag == 0)) //if left mouse button has been pressed
	{
	   gf_x_low_temp  = xflt(MouseX);
	   gf_y_high_temp = yflt(MouseY);
	   gi_zoom_flag = 1;
	}
	if ((MouseLp == TRUE) && (MouseL == TRUE) && (gi_zoom_flag == 1)) //if left mouse button not yet pressed
    {
	   frect(gf_x_low_temp,gf_y_high_temp,gf_x_move,gf_y_move,'w');
	   gf_x_move = xflt(MouseX);
	   gf_y_move = yflt(MouseY);
	   if ((gf_x_move > gf_x_low_temp) && (gf_y_move < gf_y_high_temp))
	   {
	     frect(gf_x_low_temp,gf_y_high_temp,gf_x_move,gf_y_move,'s');
	   }
	}

	if ((MouseLp == TRUE) && (MouseL == FALSE) && (gi_zoom_flag == 1)) //if left mouse button has been pressed
	{
	   gf_x_high_temp = xflt(MouseX);
	   gf_y_low_temp  = yflt(MouseY);

	   if ((gf_x_high_temp >= gf_x_low_temp) &&
           (gf_y_high_temp >= gf_y_low_temp))
	   {
		   gf_x_low  = gf_x_low_temp;
		   gf_x_high = gf_x_high_temp;
		   gf_y_low  = gf_y_low_temp;
		   gf_y_high = gf_y_high_temp;

		   Cls();//clear screen
		   axis(gf_x_low,gf_x_high,gf_y_low,gf_y_high);
	   }
	   MouseLp == FALSE;
	   gi_zoom_flag = 0;
	}

	if ((MouseRp == TRUE) && (MouseR == FALSE))//right mouse button undoes zooming
	{
	  gf_x_low  = gf_x_low_sav;
	  gf_x_high = gf_x_high_sav;
	  gf_y_low  = gf_y_low_sav;
	  gf_y_high = gf_y_high_sav;
	  Cls();//clear screen
	  axis(gf_x_low,gf_x_high,gf_y_low,gf_y_high);
      MouseRp = FALSE;
	}


	
//======Remove previous graph.============================
//------We do this simply by repainting in white.---------
        draw_graph(gfa_xmold, i_D,'w');

//======Begin with the new drawings.==================================
//------draw surrounding box------------------------------------------
		box('y');//grey box

//-----draw grid------------------------------------------------------
	    grid('y',5,5);

//-----draw tolerance scheme--------------------------------------------
		gf_ygraph = fa_bound[i_D-1];
		if (gf_y_high < gf_ygraph) gf_ygraph = gf_y_high;

		fline(-1,gf_y_high,-1,+1,'b'); 
		fline(+1,gf_y_high,+1,+1,'b');
		fline(-1,+1,+1,+1,'b'); 
		fline(-1.2,-1,+1.2,-1,'b'); 
		fline(-1.2,-1,-1.2,gf_ygraph,'b'); 
		fline(+1.2,-1,+1.2,gf_ygraph,'b'); 

//------Draw new polynomial----------------------------------------------

        draw_graph(best, i_D, 'r');

//----Save old polynomial so that it can be erased later------------------

		for (j=0; j<i_D; j++)
		{
		  gfa_xmold[j] = best[j];
		}

//------Display current number of trials----------------------------------
// x-positioning (gf_x_low + k*(gf_x_high - gf_x_low))

		 sprintf(gca_sbuft,"No. of trials:  %ld",l_nfeval);
		 myprint((gf_x_low + 0*(gf_x_high - gf_x_low)),(gf_y_low - 0.1*(gf_y_high - gf_y_low)),gca_sbuft);

//------Display current Iteration-----------------------------------------
		 sprintf(gca_sbuft,"Iteration:       %d",i_gen);
		 myprint((gf_x_low + 0.3*(gf_x_high - gf_x_low)),(gf_y_low - 0.1*(gf_y_high - gf_y_low)),gca_sbuft);

//------Display current pathlength----------------------------------------
		 sprintf(gca_sbuft,"best OFUNC value: %-11.10g",f_emin);
		 myprint((gf_x_low + 0.5*(gf_x_high - gf_x_low)),(gf_y_low - 0.1*(gf_y_high - gf_y_low)),gca_sbuft);

//------Display DE strategy-----------------------------------------------
		 sprintf(gca_sbuft,"DE-strategy No.: %d",i_strategy);
		 myprint((gf_x_low + 0.9*(gf_x_high - gf_x_low)),(gf_y_low - 0.1*(gf_y_high - gf_y_low)),gca_sbuft);

}