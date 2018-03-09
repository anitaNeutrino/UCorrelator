
void plotBaseList(){
	// int N = getNBases();
	int N = 555;
	// int* x = {0};
	// int* y = {0};
	TGraphAntarctica* plot = new TGraphAntarctica(N);
	for (int i = 0; i < N; i++){
		std::cout<<i<< " "<<  BaseList::getAbstractBase(i).getName() << std::endl;
		plot->SetPoint(i, BaseList::getAbstractBase(i).getPosition(0));
	}
	plot->Draw();
}
