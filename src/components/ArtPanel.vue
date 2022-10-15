<template>
  <v-container fluid style="max-width: 1280px">
    <v-row class="text-left">
      <v-col cols="12">
        <h1 class="secondary--text my-3">My Art</h1>
        <v-carousel
          cycle
          height="720"
          hide-delimiters
          next-icon="arrow_forward"
          prev-icon="arrow_back"
          class="art-carousel"
        >
          <v-carousel-item
            v-for="(item, i) in art"
            :key="i"
            :src="item.picture"
            contain
            style="background-color: #ffeeee"
          >
          </v-carousel-item>
        </v-carousel>
      </v-col>
    </v-row>
  </v-container>
</template>
<script>
import Vue from "vue";

export default {
  name: "ResearchPanel",
  data() {
    return {
      art: [
        {
          picture: require("../assets/aChaoticPresentation.png"),
        },
        {
          picture: require("../assets/backFromTheDead.png"),
        },
        {
          picture: require("../assets/dackel.png"),
        },
        {
          picture: require("../assets/dogAtTheEndOfTheRainbow.png"),
        },
        {
          picture: require("../assets/dontEatMe.png"),
        },
        {
          picture: require("../assets/fingersAndOcToes.png"),
        },
        {
          picture: require("../assets/forceOfNature.png"),
        },
        {
          picture: require("../assets/frosty.png"),
        },
        {
          picture: require("../assets/mightyPirate.png"),
        },
        {
          picture: require("../assets/myUnbrokenHeart.png"),
        },
        {
          picture: require("../assets/pigAtTheLavaPool.png"),
        },
        {
          picture: require("../assets/spiltMilkDay.png"),
        },
        {
          picture: require("../assets/spiltMilkNight.png"),
        },
      ],
    };
  },
  created() {
    Vue.set(this, "art", this.randomiseArrayOrder(this.art));
  },
  methods: {
    goToUrl(url) {
      window.open(url, "_blank");
    },
    randomiseArrayOrder(array) {
      let currentIndex = array.length,
        randomIndex;

      // While there remain elements to shuffle.
      while (currentIndex != 0) {
        // Pick a remaining element.
        randomIndex = Math.floor(Math.random() * currentIndex);
        currentIndex--;

        // And swap it with the current element.
        [array[currentIndex], array[randomIndex]] = [
          array[randomIndex],
          array[currentIndex],
        ];
      }

      return array;
    },
  },
};
</script>
<style lang="scss" scoped>
.art-carousel {
  height: 720px;
  border: #d93908;
  border-style: solid;
  border-width: 1px;
  border-radius: 0.25em;
  box-shadow: 0px 0px 5px #d93908;
}
</style>
