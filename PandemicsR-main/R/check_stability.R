check_stability<- function(state, t, t0){

    if (t - t0 >= 5){

      stable_now <-
        (state$opinion_changes/(t - t0)) < 0.1 &&
        (state$membership_changes/(t - t0)) < 0.1


      if(stable_now){

        state$stable_windows <-
          state$stable_windows + 1

      } else {

        state$stable_windows <- 0

      }


      state$stable <-
        state$stable_windows >= 3


      state$opinion_changes <- 0
      state$membership_changes <- 0

      t0 <- t
    }

  return(list(state, t0))
}
