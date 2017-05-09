#include "common.h"

class kinetochore
{
	protected:
		double pos;
		double prev_pos;
		double drag;
	public:
		kinetochore()
		{
			prev_pos = 0.5*LZ;
			pos = prev_pos;
			drag = KT_DRAG;
		}
		
//get variables
		double get_pos() const {return pos;}
		double get_prev_pos() const {return prev_pos;}
		double get_drag() const {return drag;}

//set variables
		void set_pos(double pos_in){pos = pos_in;}

//set prev
		void set_prev_pos()
		{
			if(pos > LZ)
			  pos = LZ;//LZ+1 --> LZ on 120323
			else if(pos < 0.)
			  pos = 0.;
			prev_pos = pos;
		}


		void load(double force_in)
		{
			pos += force_in * dt / KT_DRAG;
		}


		void move(double movement)
		{
			pos += movement;
		}		
};

extern vector<kinetochore> kt_list;

