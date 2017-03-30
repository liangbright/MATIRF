function [x y A R b A_std Region]=GaussMixFit(MethodStr, x_Init, y_Init, A_Init, R_Init, b_Init, Radius_xy, Range_A, Range_R, Range_b, ...                                                    
                                              OtherParam, I, ObjLabel, Ib_Prior, Ib_SigmaPrior, A_Prior, A_SigmaPrior)
                                                

switch MethodStr
    case 'A_b'        
        [A b A_std Region]=GaussMixFit_A_b(x_Init, y_Init, A_Init, R_Init, b_Init, Range_A, Range_b, ...                                                    
                                           OtherParam, I, ObjLabel, Ib_Prior, Ib_SigmaPrior, A_Prior, A_SigmaPrior);
        x=x_Init; y=y_Init; R=R_Init;                                     
                
    case 'x_y_A_b'
        [x y A b A_std Region] =GaussMixFit_x_y_A_b(x_Init, y_Init, A_Init, R_Init, b_Init, Radius_xy, Range_A, Range_b, ...                                                    
                                                    OtherParam, I, ObjLabel, Ib_Prior, Ib_SigmaPrior, A_Prior, A_SigmaPrior);
        R=R_Init;
                            
    case 'x_y_A_bPlane'
        [x y A b A_std Region] =GaussMixFit_x_y_A_bPlane(x_Init, y_Init, A_Init, R_Init, b_Init, Radius_xy, Range_A, Range_b, ...                                                    
                                                         OtherParam, I, ObjLabel, Ib_Prior, Ib_SigmaPrior, A_Prior, A_SigmaPrior);
        R=R_Init;
                             
    case 'x_y_A_R_b'
        [x y A R b A_std Region] =GaussMixFit_x_y_A_R_b(x_Init, y_Init, A_Init, R_Init, b_Init, Radius_xy, Range_A, Range_R, Range_b, ...                                                    
                                                        OtherParam, I, ObjLabel, Ib_Prior, Ib_SigmaPrior, A_Prior, A_SigmaPrior);
                        
    case 'x_y_A_R_bPlane'
        [x y A R b A_std Region] =GaussMixFit_x_y_A_R_bPlane(x_Init, y_Init, A_Init, R_Init, b_Init, Radius_xy, Range_A, Range_R, Range_b, ...                                                    
                                                             OtherParam, I, ObjLabel, Ib_Prior, Ib_SigmaPrior, A_Prior, A_SigmaPrior);        
                                                         
    otherwise
        disp('Wrong @ GaussMixFit')
end