<template>
  <v-container fluid>
    <v-row class="text-left">
      <v-col cols="12" class="pa-0">
        <v-carousel
          cycle
          hide-delimiters
          next-icon="arrow_forward"
          prev-icon="arrow_back"
          :height="viewPortHeight"
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
          picture: require("../assets/dumpling.png"),
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
      windowHeight: window.innerHeight,
    };
  },
  computed: {
    viewPortHeight() {
      // The tabs at the top are 48 pixels in height
      return this.windowHeight - 48;
    },
  },
  created() {
    Vue.set(this, "art", this.randomiseArrayOrder(this.art));
    window.addEventListener("resize", () => {
      this.windowHeight = window.innerHeight;
    });
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
}
</style>
