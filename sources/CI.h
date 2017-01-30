
#include<list>

void GetCI(list<double>& l, double c, double& median, double& min, double& max)	{

	l.sort();
	int n = ((int) (((double) l.size()) * (1-c)));
	list<double>::const_iterator i = l.begin();
	for (int j=0; j<n; j++)	{
		i++;
	}
	min = *i;

	n = ((int) (((double) l.size()) * c));
	i = l.begin();
	for (int j=0; j<n; j++)	{
		i++;
	}
	max = *i;

	i = l.begin();
	for (int j=0; j<((int) (l.size() / 2)); j++)	{
		i++;
	}
	median = *i;
}

