class WholePathMove{
  public:
  int col;
  long accepted;
  long total;
    real width;
    void setwidth(real wi){
      width=wi;
    }
    WholePathMove (){
      width=0.5;
      col=4;
      total=0;
      accepted=0;
      };
    int update();
    void finish();
};
class SingleSliceMove{
  public:
  int col;
  long accepted;
  long total;
    real width;
    void setwidth(real wi){
      width=wi;
    }
    SingleSliceMove(){
      width=0.5;
      col=4;
      total=0;
      accepted=0;
      };
    int update();
    void finish();
};
